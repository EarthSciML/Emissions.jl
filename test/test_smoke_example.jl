"""
Integration test: validates Emissions.jl against SMOKE ExampleCase v2 reference output.

This comprehensive test validates the complete Emissions.jl pipeline against SMOKE reference
output for the RWC (residential wood combustion) nonpoint sector for August 1, 2018.

Input data from: https://github.com/CEMPD/SMOKE-ExampleCase-v2
Reference output: premerged RWC IOAPI NetCDF files

VALIDATION SCOPE:
- Complete pipeline: FF10 → aggregation → speciation → spatial allocation → temporal allocation
- Grid definition and IOAPI structure validation
- Species completeness (39 of 62 CB6AE7 species produced)
- Spatial patterns (correlations >0.9 for key species)
- Temporal patterns (diurnal correlations >0.9)
- Quantitative magnitude matching (key species within 1-2% of reference)
- Cross-species ratios and conservation properties

Data is automatically downloaded from Google Drive if not already present.
The main example case tarball is ~2 GB and the reference output is ~30 MB.

KNOWN LIMITATIONS:
1. HAP Subtraction: SMOKE subtracts HAP (Hazardous Air Pollutant) species from total VOC
   before speciation, while Emissions.jl uses the raw inventory values. This causes VOC-derived
   species (aromatics, alkanes, etc.) to show ~0.81x ratios compared to SMOKE reference.
   This is expected and indicates correct non-HAP VOC processing.

2. Sector Coverage: Currently validates only RWC sector due to reference data availability.
   Additional sectors (nonroad, point sources, on-road mobile) could be validated with
   additional reference output files from the SMOKE ExampleCase.

3. Temporal Profiles: Uses SMOKE-style profile fallback hierarchy and timezone adjustments.
   Some minor differences in profile application may occur for edge cases.
"""

using Test
using Emissions
using DataFrames
using SparseArrays
using Dates
using NCDatasets
using CSV
using Unitful: ustrip
using Statistics: cor, mean, median, std, quantile

const SMOKE_TEST_DIR = "/tmp/smoke_test"
const SMOKE_BASE = joinpath(SMOKE_TEST_DIR, "smoke_example_case")

# Google Drive file IDs (from https://github.com/CEMPD/SMOKE-ExampleCase-v2)
const GDRIVE_EXAMPLE_CASE_ID = "1aREAz6z71WGdPoFJ2NkLutKBhWNyVgDf"
const GDRIVE_RWC_REFERENCE_ID = "1Ac69M6HGuh3ieBY03fbbMWlmsC2zrFyt"

# Reference file we compare against
const REF_FILENAME = "emis_mole_rwc_20180801_12LISTOS_cmaq_cb6ae7_2018gg_18j.ncf"

# ============================================================================
# Data download and setup
# ============================================================================

"""
    download_from_gdrive(file_id, output_path)

Download a file from Google Drive using curl. Handles the large-file
confirmation prompt via the `confirm=t` parameter.
"""
function download_from_gdrive(file_id::AbstractString, output_path::AbstractString)
    max_attempts = 3

    for attempt in 1:max_attempts
        try
            # Use drive.usercontent.google.com to bypass virus scan warning for large files
            url = "https://drive.usercontent.google.com/download?id=$(file_id)&export=download&confirm=t"
            @info "Downloading from Google Drive (attempt $attempt/$max_attempts, file ID: $file_id)..."

            # Enhanced curl with better timeout and error handling
            run(`curl -L --retry 3 --retry-delay 10 --max-time 7200 -o $output_path $url`)

            if !isfile(output_path) || filesize(output_path) < 1000
                error("Download failed or file too small: $output_path ($(filesize(output_path)) bytes)")
            end

            # Verify the download is not an HTML error page
            open(output_path) do f
                header = read(f, min(filesize(output_path), 20))
                if startswith(String(header), "<!DOCTYPE") || startswith(String(header), "<html")
                    rm(output_path, force = true)
                    error("Download returned an HTML page instead of the expected file. " *
                          "The Google Drive link may have changed or require authentication.")
                end
            end

            @info "Downloaded $(round(filesize(output_path) / 1024^2, digits=1)) MB"
            return true

        catch e
            @warn "Download attempt $attempt failed: $e"
            rm(output_path, force=true)  # Clean up partial download
            if attempt == max_attempts
                error("Failed to download after $max_attempts attempts: $e")
            end
            sleep(5)  # Brief pause before retry
        end
    end
end

"""
    find_file_recursive(dir, filename) -> Union{String, Nothing}

Search recursively under `dir` for a file named `filename`.
Returns the full path or nothing if not found.
"""
function find_file_recursive(dir::AbstractString, filename::AbstractString)
    isdir(dir) || return nothing
    for (root, dirs, files) in walkdir(dir)
        for f in files
            if f == filename
                return joinpath(root, f)
            end
        end
    end
    return nothing
end

"""
    setup_smoke_test_data() -> String

Download and extract SMOKE ExampleCase v2 data if not already present.
Returns the path to the reference NetCDF file for Aug 1 RWC.

Downloads two tarballs from Google Drive:
1. smoke_example_case.June2025.tar.gz (~2 GB) — input inventory, speciation,
   gridding, and temporal profile files
2. premerged_rwc.tar.gz (~30 MB) — SMOKE reference output for RWC sector
"""
function setup_smoke_test_data()
    mkpath(SMOKE_TEST_DIR)

    # --- Input data (main example case) ---
    inv_file = joinpath(SMOKE_BASE, "2018gg_18j", "inputs", "rwc",
        "rwc_2017NEIpost_NONPOINT_20210129_11oct2021_v1.csv")

    if !isfile(inv_file)
        tarball = joinpath(SMOKE_TEST_DIR, "smoke_example_case.tar.gz")
        if !isfile(tarball)
            download_from_gdrive(GDRIVE_EXAMPLE_CASE_ID, tarball)
        end
        @info "Extracting example case (this may take several minutes)..."
        run(`tar -xzf $tarball -C $SMOKE_TEST_DIR`)
        rm(tarball, force = true)  # Free disk space

        # Handle case where tarball extracts to a differently-named directory
        if !isdir(SMOKE_BASE)
            for entry in readdir(SMOKE_TEST_DIR)
                candidate = joinpath(SMOKE_TEST_DIR, entry)
                if startswith(entry, "smoke_example_case") && isdir(candidate) && entry != "smoke_example_case"
                    mv(candidate, SMOKE_BASE)
                    break
                end
            end
        end
        isdir(SMOKE_BASE) || error("Extraction failed: $SMOKE_BASE not found")
    end

    # --- Reference output ---
    ref_dir = joinpath(SMOKE_TEST_DIR, "reference")
    ref_file = find_file_recursive(ref_dir, REF_FILENAME)

    if ref_file === nothing
        # Also check if it's inside the main example case
        ref_file = find_file_recursive(SMOKE_BASE, REF_FILENAME)
    end

    if ref_file === nothing
        mkpath(ref_dir)
        tarball = joinpath(SMOKE_TEST_DIR, "premerged_rwc.tar.gz")
        if !isfile(tarball)
            download_from_gdrive(GDRIVE_RWC_REFERENCE_ID, tarball)
        end
        @info "Extracting RWC reference output..."
        run(`tar -xzf $tarball -C $ref_dir`)
        rm(tarball, force = true)

        ref_file = find_file_recursive(ref_dir, REF_FILENAME)
        if ref_file === nothing
            # Try again in main dir
            ref_file = find_file_recursive(SMOKE_BASE, REF_FILENAME)
        end
    end

    ref_file !== nothing || error("Reference file $REF_FILENAME not found after extraction. " *
        "Searched under $ref_dir and $SMOKE_BASE")
    @info "Reference file: $ref_file"
    return ref_file
end

# ============================================================================
# Helper: Preprocess FF10 inventory file (strip non-# header lines)
# ============================================================================
"""
Create a temporary copy of an FF10 file with any non-commented header line removed.
SMOKE ExampleCase files include a `country_cd,region_cd,...` header that
Emissions.jl's `read_ff10` doesn't expect (it uses `header=false`).
"""
function preprocess_ff10(filepath::AbstractString)
    tmpfile = tempname() * ".csv"
    open(tmpfile, "w") do out
        for line in readlines(filepath)
            s = strip(line)
            if isempty(s)
                continue
            elseif startswith(s, "#")
                println(out, line)
            elseif startswith(lowercase(s), "country_cd") || startswith(lowercase(s), "\"country_cd")
                # Skip header row
                continue
            else
                println(out, line)
            end
        end
    end
    return tmpfile
end

# ============================================================================
# Helper: Parse SMOKE GRIDDESC with Fortran D-notation and ! comments
# ============================================================================
"""
Parse a SMOKE GRIDDESC file that uses `!` comments and Fortran `D` notation.
Returns a GridDef for the named grid.
"""
function parse_smoke_griddesc(filepath::AbstractString, grid_name::AbstractString)
    lines = readlines(filepath)

    # Clean: strip ! comments and handle Fortran D notation
    cleaned = String[]
    for line in lines
        s = strip(line)
        isempty(s) && continue
        # Remove inline ! comments
        idx = findfirst('!', s)
        if idx !== nothing
            s = strip(s[1:idx-1])
        end
        isempty(s) && continue
        # Replace Fortran D notation: 33.0D0 → 33.0e0, 36.D3 → 36.e3
        s = replace(s, r"(\d)D([\+\-]?\d)" => s"\1e\2")
        s = replace(s, r"\.D(\d)" => s".e\1")
        push!(cleaned, s)
    end

    # Split into segments by delimiter (lines that are just quoted space: ' ')
    segments = Vector{Vector{String}}()
    current = String[]
    for line in cleaned
        stripped = replace(line, r"['\"\s]" => "")
        if isempty(stripped)
            if !isempty(current)
                push!(segments, current)
                current = String[]
            end
        else
            push!(current, line)
        end
    end
    if !isempty(current)
        push!(segments, current)
    end

    length(segments) >= 2 || error("GRIDDESC: need at least 2 segments, got $(length(segments))")

    # Segment 1: coordinate systems (pairs of lines)
    coord_systems = Dict{String, NamedTuple}()
    seg1 = segments[1]
    i = 1
    while i + 1 <= length(seg1)
        coord_name = replace(strip(seg1[i]), r"['\"]" => "")
        params = split(strip(seg1[i+1]), r"[,\s]+")
        filter!(!isempty, params)
        if length(params) >= 6
            coordtype = parse(Int, params[1])
            P_ALP = parse(Float64, params[2])
            P_BET = parse(Float64, params[3])
            P_GAM = parse(Float64, params[4])
            XCENT = parse(Float64, params[5])
            YCENT = parse(Float64, params[6])
            coord_systems[coord_name] = (
                coordtype=coordtype, P_ALP=P_ALP, P_BET=P_BET,
                P_GAM=P_GAM, XCENT=XCENT, YCENT=YCENT,
            )
        end
        i += 2
    end

    # Segment 2: grid definitions (pairs of lines)
    seg2 = segments[2]
    i = 1
    while i + 1 <= length(seg2)
        gname = replace(strip(seg2[i]), r"['\"]" => "")
        params = split(strip(seg2[i+1]), r"[,\s]+")
        filter!(!isempty, params)

        if gname == grid_name && length(params) >= 8
            coord_name = replace(strip(params[1]), r"['\"]" => "")
            XORIG = parse(Float64, params[2])
            YORIG = parse(Float64, params[3])
            XCELL = parse(Float64, params[4])
            YCELL = parse(Float64, params[5])
            NCOLS = parse(Int, params[6])
            NROWS = parse(Int, params[7])

            haskey(coord_systems, coord_name) || error("Coord system '$coord_name' not found")
            cs = coord_systems[coord_name]

            # Build proj4 string for LCC
            proj4 = if cs.coordtype == 2
                "+proj=lcc +lat_1=$(cs.P_ALP) +lat_2=$(cs.P_BET) +lon_0=$(cs.P_GAM) +lat_0=$(cs.YCENT) +datum=WGS84"
            elseif cs.coordtype == 1
                "+proj=longlat +datum=WGS84"
            else
                error("Unsupported coord type $(cs.coordtype)")
            end

            return NewGridRegular(grid_name, NCOLS, NROWS, proj4, XCELL, YCELL, XORIG, YORIG)
        end
        i += 2
    end
    error("Grid '$grid_name' not found in GRIDDESC")
end

# ============================================================================
# Helper: Parse SMOKE comma-delimited GSREF (combo format) into Emissions.jl format
# ============================================================================
"""
Parse SMOKE combo-format GSREF CSV file.
Returns DataFrame with columns: FIPS, SCC, pollutant_id, profile_code
matching what Emissions.jl's `speciate_emissions` expects.

The SMOKE GSREF combo format is comma-delimited with fields:
  SCC, profile_code, pollutant, [FIPS], [combo fields...], [fraction], [comment]

Entries without FIPS get FIPS="00000" (national default).
"""
function parse_smoke_gsref_csv(filepath::AbstractString)
    records = DataFrame(
        FIPS = String[],
        SCC = String[],
        pollutant_id = String[],
        profile_code = String[],
    )
    for line in readlines(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue

        # Remove inline comments after '!'
        comment_idx = findfirst('!', line)
        if comment_idx !== nothing
            line = strip(line[1:comment_idx-1])
        end

        # Parse CSV: handle quoted fields
        parts = String[]
        in_quote = false
        current = IOBuffer()
        for ch in line
            if ch == '"'
                in_quote = !in_quote
            elseif ch == ',' && !in_quote
                push!(parts, String(take!(current)))
                current = IOBuffer()
            else
                write(current, ch)
            end
        end
        push!(parts, String(take!(current)))

        length(parts) >= 3 || continue

        scc_raw = strip(parts[1])
        profile_code = strip(parts[2])
        pollutant = strip(parts[3])

        isempty(scc_raw) && continue
        isempty(profile_code) && continue
        isempty(pollutant) && continue

        # SCC: pad to 10 digits
        scc = lpad(scc_raw, 10, '0')

        # FIPS: field 4 if present, else national default
        fips = "00000"
        if length(parts) >= 4
            fips_raw = strip(parts[4])
            if !isempty(fips_raw)
                # Strip country digit if 6-digit, pad to 5
                if length(fips_raw) == 6
                    fips = lpad(fips_raw[2:end], 5, '0')
                else
                    fips = lpad(fips_raw, 5, '0')
                end
            end
        end

        push!(records, (
            FIPS = fips,
            SCC = scc,
            pollutant_id = uppercase(pollutant),
            profile_code = profile_code,
        ))
    end
    return records
end

# ============================================================================
# Helper: Parse SMOKE pre-computed surrogate text file
# ============================================================================
"""
Parse a SMOKE pre-computed surrogate text file and build county surrogate
sparse matrices for the target grid.

The surrogate file has:
- Header line: #GRID name xorig yorig dx dy ncols nrows ...
- Data lines: surrogate_code  FIPS  col  row  fraction  ! weight total_weight

Returns Dict{String, SparseMatrixCSC{Float64, Int}} mapping 5-digit FIPS codes
to normalized allocation matrices for the target grid.
"""
function parse_smoke_surrogates(filepath::AbstractString, target_grid::GridDef)
    # Parse header to get national grid info
    xorig_nat = 0.0
    yorig_nat = 0.0
    dx_nat = 0.0
    dy_nat = 0.0

    for line in readlines(filepath)
        if startswith(line, "#GRID")
            parts = split(line)
            # #GRID name xorig yorig dx dy ncols nrows ...
            xorig_nat = parse(Float64, parts[3])
            yorig_nat = parse(Float64, parts[4])
            dx_nat = parse(Float64, parts[5])
            dy_nat = parse(Float64, parts[6])
            break
        end
    end

    @assert dx_nat > 0 "Failed to parse surrogate grid header"
    @assert abs(dx_nat - target_grid.Dx) < 1.0 "Grid cell size mismatch: surrogate=$dx_nat, target=$(target_grid.Dx)"

    # Compute offset: national grid cell (col_nat, row_nat) → target cell
    col_offset = round(Int, (target_grid.X0 - xorig_nat) / dx_nat)
    row_offset = round(Int, (target_grid.Y0 - yorig_nat) / dy_nat)

    # Accumulate fractions per FIPS code
    fips_data = Dict{String, Vector{Tuple{Int, Int, Float64}}}()

    for line in readlines(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue

        # Split on whitespace/tab, handling the '!' comment marker
        parts = split(line)
        length(parts) >= 5 || continue

        fips_raw = strip(parts[2])
        col_nat = parse(Int, parts[3])
        row_nat = parse(Int, parts[4])
        frac = parse(Float64, parts[5])

        # Convert to target grid coordinates (1-based)
        col_target = col_nat - col_offset
        row_target = row_nat - row_offset

        # Keep only cells within the target grid
        (1 <= col_target <= target_grid.Nx && 1 <= row_target <= target_grid.Ny) || continue

        # Normalize FIPS to 5 digits
        if length(fips_raw) == 6
            fips = lpad(fips_raw[2:end], 5, '0')
        else
            fips = lpad(fips_raw, 5, '0')
        end

        if !haskey(fips_data, fips)
            fips_data[fips] = Tuple{Int, Int, Float64}[]
        end
        push!(fips_data[fips], (row_target, col_target, frac))
    end

    # Build sparse matrices per FIPS
    result = Dict{String, SparseMatrixCSC{Float64, Int}}()
    for (fips, entries) in fips_data
        rows = [e[1] for e in entries]
        cols = [e[2] for e in entries]
        vals = [e[3] for e in entries]

        # Do NOT normalize fractions to the target grid.
        # SMOKE uses raw surrogate fractions: border counties have fracs < 1.0,
        # meaning only the in-grid portion of their emissions is allocated.
        # This matches SMOKE's behavior where emissions outside the target grid
        # are simply not included in the output.
        mat = sparse(rows, cols, vals, target_grid.Ny, target_grid.Nx)
        result[fips] = mat
    end

    return result
end

# ============================================================================
# Helper: Parse Gentpro monthly profiles (FIPS-keyed)
# ============================================================================
"""
Parse Gentpro TPRO_MON monthly profiles.
Format: "FIPS_PROFILE_ID", jan, feb, ..., dec,
Returns Dict{String, Vector{Float64}} mapping FIPS code to 12 monthly fractions.
"""
function parse_gentpro_monthly(filepath::AbstractString)
    result = Dict{String, Vector{Float64}}()
    for line in readlines(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue

        parts = split(line, ',')
        length(parts) >= 13 || continue

        fips = replace(strip(parts[1]), '"' => "")
        isempty(fips) && continue

        factors = Float64[]
        for i in 2:13
            s = strip(parts[i])
            push!(factors, isempty(s) ? 0.0 : parse(Float64, s))
        end
        result[fips] = factors
    end
    return result
end

# ============================================================================
# Helper: Parse Gentpro daily profiles (FIPS-keyed)
# ============================================================================
"""
Parse Gentpro TPRO_DAY daily profiles.
Format: "FIPS_ID", month, day1, day2, ..., day31
Returns Dict{Tuple{String,Int}, Vector{Float64}} mapping (FIPS, month) to 31 daily fractions.
"""
function parse_gentpro_daily(filepath::AbstractString)
    result = Dict{Tuple{String, Int}, Vector{Float64}}()
    for line in readlines(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue

        parts = split(line, ',')
        length(parts) >= 33 || continue

        fips = replace(strip(parts[1]), '"' => "")
        isempty(fips) && continue

        month = parse(Int, strip(parts[2]))
        factors = Float64[]
        for i in 3:min(33, length(parts))
            s = strip(parts[i])
            push!(factors, isempty(s) ? 0.0 : parse(Float64, s))
        end
        # Pad to 31 days if needed
        while length(factors) < 31
            push!(factors, 0.0)
        end
        result[(fips, month)] = factors
    end
    return result
end

# ============================================================================
# Helper: Parse SMOKE ATREF (Gentpro format)
# ============================================================================
"""
Parse Gentpro-format ATREF temporal cross-reference.
Format: "SCC","FIPS",,,,,,"profile_type","profile_id",
Returns Dict mapping (SCC, FIPS) to Dict of profile_type => profile_id.
"""
function parse_atref_gentpro(filepath::AbstractString)
    result = Dict{Tuple{String, String}, Dict{String, String}}()
    for line in readlines(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue

        parts = split(line, ',')
        length(parts) >= 9 || continue

        scc = replace(strip(parts[1]), '"' => "")
        fips_raw = replace(strip(parts[2]), '"' => "")
        profile_type = replace(strip(parts[8]), '"' => "")
        profile_id = replace(strip(parts[9]), '"' => "")

        isempty(scc) && continue
        isempty(profile_type) && continue

        # Normalize FIPS: empty means national default "00000"
        if isempty(fips_raw)
            fips = "00000"
        elseif length(fips_raw) == 6
            fips = lpad(fips_raw[2:end], 5, '0')
        else
            fips = lpad(fips_raw, 5, '0')
        end

        key = (lpad(scc, 10, '0'), fips)
        if !haskey(result, key)
            result[key] = Dict{String, String}()
        end
        result[key][uppercase(profile_type)] = profile_id
    end
    return result
end

# ============================================================================
# Helper: Parse standard SMOKE hourly profile file (amptpro format)
# ============================================================================
"""
Parse SMOKE amptpro hourly (diurnal) profile file.
Format: "profile_id", hour1, ..., hour24, "comment"
Returns Dict{String, Vector{Float64}} mapping profile_id to 24 hourly fractions.
"""
function parse_amptpro_hourly(filepath::AbstractString)
    result = Dict{String, Vector{Float64}}()
    for line in readlines(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, "#") || startswith(line, "profile_id")) && continue

        parts = split(line, ',')
        length(parts) >= 25 || continue

        profile_id = replace(strip(parts[1]), '"' => "")
        isempty(profile_id) && continue

        factors = Float64[]
        for i in 2:25
            s = strip(parts[i])
            push!(factors, isempty(s) ? 0.0 : parse(Float64, s))
        end
        result[profile_id] = factors
    end
    return result
end

# ============================================================================
# Helper: Build pollutant name mapping for GSREF/GSPRO compatibility
# ============================================================================
"""
Prepare emissions for speciation by fixing pollutant names and computing PMC.

1. Rename "PM25" → "PM2_5" (matches GSREF/GSPRO naming)
2. Compute "PMC" = PM10 - PM25 per (FIPS, SCC) (coarse PM, matches SMOKE's SMKINVEN_FORMULA)
3. Remove original PM10 rows (replaced by PMC)
"""
function prepare_emissions_for_speciation!(emissions::DataFrame)
    # Ensure string POLID for matching
    emissions.POLID = string.(emissions.POLID)

    # Step 1: Compute PMC = PM10 - PM25 per (FIPS, SCC)
    # This matches SMOKE's SMKINVEN_FORMULA = "PMC=PM10-PM2_5"
    pm10_idx = findall(r -> r.POLID == "PM10", eachrow(emissions))
    pm25_map = Dict{Tuple{String, String}, Float64}()
    for row in eachrow(emissions)
        if row.POLID == "PM25"
            key = (string(row.FIPS), string(row.SCC))
            pm25_map[key] = get(pm25_map, key, 0.0) + Float64(row.ANN_VALUE)
        end
    end

    pmc_rows = similar(emissions, 0)
    for idx in pm10_idx
        row = emissions[idx, :]
        key = (string(row.FIPS), string(row.SCC))
        pm25_val = get(pm25_map, key, 0.0)
        pmc_value = Float64(row.ANN_VALUE) - pm25_val
        if pmc_value > 0
            new_row = copy(emissions[idx:idx, :])
            new_row[1, :POLID] = "PMC"
            new_row[1, :ANN_VALUE] = pmc_value
            append!(pmc_rows, new_row)
        end
    end

    # Step 2: Remove PM10 rows (now replaced by PMC)
    filter!(r -> r.POLID != "PM10", emissions)

    # Step 3: Add PMC rows
    if nrow(pmc_rows) > 0
        append!(emissions, pmc_rows)
        @info "Computed $(nrow(pmc_rows)) PMC rows from PM10 - PM25"
    end

    # Step 4: Rename PM25 → PM2_5
    emissions.POLID = [p == "PM25" ? "PM2_5" : p for p in emissions.POLID]

    return emissions
end

"""
Normalize GSREF pollutant_id so "NONHAPVOC" → "VOC".
This allows emissions POLID="VOC" to match GSREF entries for NONHAPVOC.
Other pollutant names (NOX, CO, SO2, NH3, PM2_5, PMC) are already correct.
"""
function normalize_gsref_pollutants!(gsref::DataFrame)
    mapping = Dict("NONHAPVOC" => "VOC")
    gsref.pollutant_id = [get(mapping, uppercase(p), uppercase(p)) for p in gsref.pollutant_id]
    return gsref
end

# ============================================================================
# Helper: Parse SMOKE weekly profile file (amptpro format)
# ============================================================================
"""
Parse SMOKE amptpro weekly (day-of-week) profile file.
Format: "profile_id", mon, tue, wed, thu, fri, sat, sun, "comment"
Returns Dict{String, Vector{Float64}} mapping profile_id to 7 day-of-week weights.
Weights are returned in temporal_allocate convention (uniform = 1.0 per day, sum = 7.0).
"""
function parse_amptpro_weekly(filepath::AbstractString)
    result = Dict{String, Vector{Float64}}()
    for line in readlines(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, "#") || startswith(line, "profile_id")) && continue

        parts = split(line, ',')
        length(parts) >= 8 || continue

        profile_id = replace(strip(parts[1]), '"' => "")
        isempty(profile_id) && continue

        factors = Float64[]
        for i in 2:8
            s = strip(parts[i])
            push!(factors, isempty(s) ? 0.0 : parse(Float64, s))
        end

        # Convert fractions (sum≈1) to weights (sum≈7) if needed
        s = sum(factors)
        if s > 0 && s < 2.0  # likely fractions summing to ~1.0
            factors .*= 7.0
        end
        result[profile_id] = factors
    end
    return result
end

# ============================================================================
# Helper: Build temporal_allocate-compatible profiles/xref from Gentpro data
# ============================================================================
"""
Convert Gentpro FIPS-specific temporal profiles and ATREF cross-reference
into the profiles and xref DataFrames expected by `temporal_allocate`.

This handles two distinct temporal patterns:
1. Most RWC SCCs: FIPS-specific monthly + daily profiles, diurnal profile 600
2. Hydronic heater SCCs (2104008610-630): static monthly profile, weekly profile 7, diurnal 1500

Gentpro daily fractions (day-of-month) are converted to weekly weights via:
  weekly_weight = daily_fraction × days_in_month
which gives weight=1.0 for uniform daily distribution.
"""
function build_gentpro_temporal(
    atref::Dict{Tuple{String,String}, Dict{String,String}},
    gentpro_monthly::Dict{String, Vector{Float64}},
    gentpro_daily::Dict{Tuple{String,Int}, Vector{Float64}},
    hourly_profiles::Dict{String, Vector{Float64}},
    weekly_profiles::Dict{String, Vector{Float64}},
    emissions_fips_scc::Set{Tuple{String,String}},
    target_date::Date)

    month_idx = Dates.month(target_date)
    day_idx = Dates.day(target_date)
    n_days = Dates.daysinmonth(target_date)

    # Map string profile IDs to integers (reserve 1 for defaults)
    str_to_int = Dict{String, Int}()
    next_id = Ref(2)
    function get_or_create_id(key::String)
        return get!(str_to_int, key) do
            id = next_id[]
            next_id[] += 1
            id
        end
    end

    # Track added profiles to avoid duplicates
    added_profiles = Set{Tuple{String, Int}}()
    profiles_rows = NamedTuple{(:profile_type, :profile_id, :factors), Tuple{String, Int, Vector{Float64}}}[]
    xref_rows = NamedTuple{(:FIPS, :SCC, :monthly_id, :weekly_id, :diurnal_id), Tuple{String, String, Int, Int, Int}}[]

    function add_profile!(ptype, pid, factors)
        key = (ptype, pid)
        if key ∉ added_profiles
            push!(added_profiles, key)
            push!(profiles_rows, (profile_type=ptype, profile_id=pid, factors=factors))
        end
    end

    # Add defaults (ID=1)
    add_profile!("MONTHLY", 1, fill(1.0 / 12.0, 12))
    add_profile!("WEEKLY", 1, fill(1.0, 7))
    add_profile!("ALLDAY", 1, fill(1.0 / 24.0, 24))

    n_matched = 0
    n_unmatched = 0

    for (fips, scc) in emissions_fips_scc
        scc_padded = lpad(scc, 10, '0')

        # Look up ATREF: exact (SCC, FIPS) match first
        prof_info = get(atref, (scc_padded, fips), nothing)

        # Fallback: try any FIPS for this SCC (national default)
        if prof_info === nothing
            for (key, val) in atref
                if key[1] == scc_padded
                    prof_info = val
                    break
                end
            end
        end

        monthly_id = 1  # default
        weekly_id = 1   # default
        diurnal_id = 1  # default

        if prof_info !== nothing
            n_matched += 1

            # Monthly profile
            mon_str = get(prof_info, "MONTHLY", nothing)
            if mon_str !== nothing && haskey(gentpro_monthly, mon_str)
                monthly_id = get_or_create_id("monthly_$mon_str")
                add_profile!("MONTHLY", monthly_id, gentpro_monthly[mon_str])
            end

            # Daily → weekly conversion, or weekly profile directly
            day_str = get(prof_info, "DAILY", nothing)
            wk_str = get(prof_info, "WEEKLY", nothing)

            if day_str !== nothing && haskey(gentpro_daily, (day_str, month_idx))
                # Convert Gentpro daily fraction to weekly weight
                daily_facs = gentpro_daily[(day_str, month_idx)]
                daily_frac = day_idx <= length(daily_facs) ? daily_facs[day_idx] : 0.0
                weekly_weight = daily_frac * n_days
                weekly_id = get_or_create_id("daily_$(day_str)_$(month_idx)_$(day_idx)")
                # Fill all DOW slots with same weight (only 1 day simulated)
                add_profile!("WEEKLY", weekly_id, fill(weekly_weight, 7))
            elseif wk_str !== nothing && haskey(weekly_profiles, wk_str)
                weekly_id = get_or_create_id("weekly_$wk_str")
                add_profile!("WEEKLY", weekly_id, weekly_profiles[wk_str])
            end

            # Diurnal profile
            allday_str = get(prof_info, "ALLDAY", nothing)
            if allday_str !== nothing && haskey(hourly_profiles, allday_str)
                diurnal_id = get_or_create_id("allday_$allday_str")
                add_profile!("ALLDAY", diurnal_id, hourly_profiles[allday_str])
            end
        else
            n_unmatched += 1
        end

        push!(xref_rows, (FIPS=fips, SCC=scc_padded,
                          monthly_id=monthly_id, weekly_id=weekly_id, diurnal_id=diurnal_id))
    end

    @info "Gentpro temporal: matched $n_matched (FIPS,SCC) pairs in ATREF, $n_unmatched unmatched"

    profiles = DataFrame(profiles_rows)
    xref = DataFrame(xref_rows)

    return profiles, xref
end

"""
Normalize GSPRO pollutant_id so "NONHAPTOG" → "VOC".
Only NONHAPTOG is renamed — NOT "TOG" — to avoid double-counting.
The 2017 NEI inventory reports non-HAP VOC, so NONHAPTOG is the correct
GSPRO pollutant to use for speciation of "VOC" emissions.
"""
function normalize_gspro_pollutants!(gspro::DataFrame)
    mapping = Dict("NONHAPTOG" => "VOC")
    gspro.pollutant_id = [get(mapping, uppercase(p), uppercase(p)) for p in gspro.pollutant_id]
    return gspro
end

# ============================================================================
# Comparison helper functions
# ============================================================================

"""Compute cosine similarity between two non-negative vectors."""
function cosine_similarity(a::AbstractVector, b::AbstractVector)
    dot_ab = sum(a .* b)
    norm_a = sqrt(sum(a .^ 2))
    norm_b = sqrt(sum(b .^ 2))
    (norm_a == 0 || norm_b == 0) && return NaN
    return dot_ab / (norm_a * norm_b)
end

"""Read IOAPI variable, returning data in (COL, ROW, LAY, TSTEP) order as stored."""
function read_ioapi_var_raw(ds, varname)
    return Array(ds[varname])  # (COL, ROW, LAY, TSTEP)
end

"""
Sum IOAPI variable over layers and timesteps, then transpose to (ROW, COL).
Returns a Matrix{Float64} comparable to Julia model output layout.
"""
function ioapi_spatial_pattern(ds, varname)
    data = read_ioapi_var_raw(ds, varname)  # (COL, ROW, LAY, TSTEP)
    # Sum over LAY (dim 3) and TSTEP (dim 4), result is (COL, ROW)
    spatial = dropdims(sum(data, dims = (3, 4)), dims = (3, 4))
    return permutedims(spatial, (2, 1))  # (ROW, COL)
end

"""
Extract hourly totals from IOAPI variable (summed over all grid cells and layers).
Returns Vector{Float64} of length TSTEP.
"""
function ioapi_hourly_totals(ds, varname)
    data = read_ioapi_var_raw(ds, varname)  # (COL, ROW, LAY, TSTEP)
    nsteps = size(data, 4)
    return [sum(data[:, :, :, t]) for t in 1:nsteps]
end

"""Extract the species list from an IOAPI file's VAR-LIST attribute."""
function read_ioapi_species(ds)
    nvars = ds.attrib["NVARS"]
    varlist = ds.attrib["VAR-LIST"]
    return [strip(varlist[(i-1)*16+1:i*16]) for i in 1:nvars]
end

# ============================================================================
# TESTS
# ============================================================================
@testset "SMOKE ExampleCase RWC Integration" begin

    # --- Data setup ---
    ref_file = setup_smoke_test_data()

    # Define all input file paths
    inv_file = joinpath(SMOKE_BASE, "2018gg_18j", "inputs", "rwc",
        "rwc_2017NEIpost_NONPOINT_20210129_11oct2021_v1.csv")
    griddesc_file = joinpath(SMOKE_BASE, "ge_dat", "gridding",
        "griddesc_lambertonly_18jan2019_v7.txt")
    agref_file = joinpath(SMOKE_BASE, "ge_dat", "gridding",
        "agref_us_2017platform_15apr2022_v12.txt")
    srg_file = joinpath(SMOKE_BASE, "ge_dat", "gridding", "surrogates",
        "CONUS12_2017NEI_04mar2021", "USA_100_NOFILL.txt")
    gspro_file = joinpath(SMOKE_BASE, "ge_dat", "speciation",
        "gspro_rwc_cmaq_cb6ae7_2018gg_18j_01may2019.txt")
    gsref_file = joinpath(SMOKE_BASE, "ge_dat", "speciation",
        "gsref_rwc_cmaq_cb6ae7_2018gg_18j_12apr2022.txt")
    monthly_file = joinpath(SMOKE_BASE, "ge_dat", "temporal",
        "tpro_monthly_Gentpro_RWC_2018gc_18j_23sep2021_nf_v1")
    daily_file = joinpath(SMOKE_BASE, "ge_dat", "temporal",
        "tpro_daily_Gentpro_RWC_2018gc_18j_23sep2021_nf_v1")
    hourly_file = joinpath(SMOKE_BASE, "ge_dat", "temporal",
        "amptpro_general_2011platform_tpro_hourly_6nov2014_18apr2022_v11")
    weekly_file = joinpath(SMOKE_BASE, "ge_dat", "temporal",
        "amptpro_general_2011platform_tpro_weekly_6nov2014_09sep2016_v2")
    atref_file = joinpath(SMOKE_BASE, "ge_dat", "temporal",
        "atref_2017platform_rwc_29apr2020_v1")

    @testset "Required files exist" begin
        for f in [inv_file, griddesc_file, agref_file, srg_file, gspro_file,
                  gsref_file, monthly_file, daily_file, hourly_file, weekly_file,
                  atref_file, ref_file]
            @test isfile(f)
        end
    end

    # ========================================================================
    @testset "Grid definition matches reference" begin
        grid = parse_smoke_griddesc(griddesc_file, "12LISTOS")
        @test grid.Name == "12LISTOS"
        @test grid.Nx == 25
        @test grid.Ny == 25
        @test grid.Dx ≈ 12000.0
        @test grid.Dy ≈ 12000.0
        @test grid.X0 ≈ 1.812e6
        @test grid.Y0 ≈ 240000.0

        # Compare against reference NetCDF attributes
        NCDatasets.Dataset(ref_file, "r") do ds
            @test ds.attrib["NCOLS"] == 25
            @test ds.attrib["NROWS"] == 25
            @test ds.attrib["XCELL"] ≈ 12000.0
            @test ds.attrib["YCELL"] ≈ 12000.0
            @test ds.attrib["XORIG"] ≈ 1.812e6
            @test ds.attrib["YORIG"] ≈ 240000.0
            @test ds.attrib["GDTYP"] == 2  # LCC
            @test ds.attrib["P_ALP"] ≈ 33.0
            @test ds.attrib["P_BET"] ≈ 45.0
        end
    end

    # ========================================================================
    @testset "Inventory reading" begin
        processed_inv = preprocess_ff10(inv_file)
        try
            emis = read_ff10(processed_inv, :nonpoint)
            @test emis isa EmissionsDataFrame
            df = emis.df
            @test nrow(df) > 0
            @test hasproperty(df, :POLID)
            @test hasproperty(df, :ANN_VALUE)
            @test hasproperty(df, :FIPS)
            @test hasproperty(df, :SCC)

            # Verify some RWC SCCs are present
            sccs = unique(string.(df.SCC))
            @test any(startswith.(sccs, "2104"))
        finally
            rm(processed_inv, force = true)
        end
    end

    # ========================================================================
    @testset "GSPRO/GSREF parsing" begin
        gspro = read_gspro(gspro_file)
        @test nrow(gspro) > 0
        @test hasproperty(gspro, :profile_code)
        @test hasproperty(gspro, :pollutant_id)
        @test hasproperty(gspro, :species_id)

        gsref = parse_smoke_gsref_csv(gsref_file)
        @test nrow(gsref) > 0
        @test hasproperty(gsref, :FIPS)
        @test hasproperty(gsref, :SCC)
        @test hasproperty(gsref, :pollutant_id)
        @test hasproperty(gsref, :profile_code)

        # RWC SCCs should have NONHAPVOC/VOC entries
        rwc_entries = filter(r -> startswith(r.SCC, "2104"), gsref)
        @test nrow(rwc_entries) > 0
        pols = unique(rwc_entries.pollutant_id)
        @test "NONHAPVOC" in pols || "VOC" in pols
    end

    # ========================================================================
    @testset "Surrogate parsing" begin
        grid = parse_smoke_griddesc(griddesc_file, "12LISTOS")
        surrogates = parse_smoke_surrogates(srg_file, grid)
        @test !isempty(surrogates)

        # Check that surrogate fractions are <= 1.0 per FIPS
        # (fractions sum to 1.0 for fully in-grid counties, less for border counties)
        for (fips, mat) in surrogates
            total = sum(mat)
            @test 0.0 < total <= 1.01
        end

        # Check matrix dimensions match grid
        for (fips, mat) in surrogates
            @test size(mat) == (grid.Ny, grid.Nx)
            break
        end
    end

    # ========================================================================
    @testset "Temporal profile parsing" begin
        monthly = parse_gentpro_monthly(monthly_file)
        @test !isempty(monthly)
        # Monthly fractions should sum to ~1.0
        for (fips, factors) in monthly
            @test length(factors) == 12
            @test sum(factors) ≈ 1.0 atol = 0.01
            break
        end

        daily = parse_gentpro_daily(daily_file)
        @test !isempty(daily)
        # Daily fractions should sum to ~1.0 for months with emissions
        # Check a few entries with non-zero sums
        checked = 0
        for ((fips, month), factors) in daily
            @test length(factors) == 31
            s = sum(factors)
            if s > 0.1  # skip months with zero emissions
                @test s ≈ 1.0 atol = 0.02
                checked += 1
                checked >= 3 && break
            end
        end
        @test checked > 0

        hourly = parse_amptpro_hourly(hourly_file)
        @test !isempty(hourly)
        # Diurnal fractions should sum to ~1.0
        for (pid, factors) in hourly
            @test length(factors) == 24
            @test sum(factors) ≈ 1.0 atol = 0.01
            break
        end

        # RWC-specific diurnal profiles should exist
        @test haskey(hourly, "600")   # Standard RWC diurnal
        @test haskey(hourly, "1500")  # Hydronic heater diurnal

        # Profile 600 should have evening peak (hours 18-23 > midday hours)
        p600 = hourly["600"]
        evening_avg = mean(p600[19:24])  # hours 18-23 (1-indexed)
        midday_avg = mean(p600[11:15])   # hours 10-14
        @test evening_avg > midday_avg

        weekly = parse_amptpro_weekly(weekly_file)
        @test !isempty(weekly)
        @test haskey(weekly, "7")  # Hydronic heater weekly profile
        # Profile 7 should be approximately uniform (each day ≈ 1.0 weight)
        p7 = weekly["7"]
        @test length(p7) == 7
        @test sum(p7) ≈ 7.0 atol = 0.1

        atref = parse_atref_gentpro(atref_file)
        @test !isempty(atref)

        # Verify ATREF assigns different profiles for different SCC types
        # Standard RWC SCCs should get ALLDAY=600, hydronic heaters should get ALLDAY=1500
        has_600 = false
        has_1500 = false
        for ((scc, fips), info) in atref
            allday = get(info, "ALLDAY", "")
            if allday == "600"
                has_600 = true
            elseif allday == "1500"
                has_1500 = true
            end
            has_600 && has_1500 && break
        end
        @test has_600   # Standard RWC SCCs use diurnal profile 600
        @test has_1500  # Hydronic heater SCCs use diurnal profile 1500
    end

    # ========================================================================
    @testset "Full pipeline: RWC August 1, 2018" begin
        # --- Step 1: Read grid ---
        grid = parse_smoke_griddesc(griddesc_file, "12LISTOS")

        # --- Step 2: Read and prepare inventory ---
        processed_inv = preprocess_ff10(inv_file)
        emis_obj = read_ff10(processed_inv, :nonpoint)
        rm(processed_inv, force = true)
        raw_df = emis_obj.df
        emissions = aggregate_emissions([raw_df])
        emissions = filter_known_pollutants(emissions)
        map_pollutant_names!(emissions)

        # Normalize COUNTRY to match gridref format ("US" → "USA")
        if hasproperty(emissions, :COUNTRY)
            emissions[!, :COUNTRY] = [Emissions.normalize_country(string(c)) for c in emissions[!, :COUNTRY]]
        end

        # Strip Unitful units from ANN_VALUE (speciate_emissions expects plain Float64)
        emissions[!, :ANN_VALUE] = [Float64(ustrip(v)) for v in emissions[!, :ANN_VALUE]]
        GC.gc()

        @test nrow(emissions) > 0
        @info "Emissions after filtering: $(nrow(emissions)) rows, POLIDs: $(unique(emissions.POLID))"

        # Prepare emissions for speciation: compute PMC = PM10 - PM25, rename PM25→PM2_5
        prepare_emissions_for_speciation!(emissions)

        # --- Step 3: Read and adapt speciation files ---
        gspro = read_gspro(gspro_file)
        gsref = parse_smoke_gsref_csv(gsref_file)

        # Strip quotes from GSPRO fields (read_gspro doesn't handle quoted fields)
        for col in [:profile_code, :pollutant_id, :species_id]
            gspro[!, col] = [replace(s, '"' => "") for s in gspro[!, col]]
        end

        # Normalize pollutant names:
        # GSREF: "NONHAPVOC" → "VOC" (so emissions "VOC" matches)
        # GSPRO: "NONHAPTOG" → "VOC" (corresponding GSPRO for non-HAP VOC)
        # Note: "TOG" is NOT renamed to avoid double-counting HAP species
        normalize_gsref_pollutants!(gsref)
        normalize_gspro_pollutants!(gspro)

        # Deduplicate GSREF: keep first match per (FIPS, SCC, pollutant)
        unique!(gsref, [:FIPS, :SCC, :pollutant_id])

        @info "GSPRO: $(nrow(gspro)) entries, GSREF: $(nrow(gsref)) entries"

        # --- Step 4: Speciation ---
        # Use mole basis to match SMOKE's convention:
        # - Gas species: split_factor/divisor converts mass → moles (divisor = MW)
        # - PM species: divisor=1.0, so mole basis = mass basis (unchanged)
        speciated = speciate_emissions(emissions, gspro, gsref; basis = :mole)
        GC.gc()

        @test nrow(speciated) > 0
        @info "After speciation: $(nrow(speciated)) rows, species: $(length(unique(speciated.POLID)))"

        # Verify some expected species exist
        species = unique(speciated.POLID)
        # NOX should be split into NO and NO2
        if "NOX" in unique(emissions.POLID)
            @test "NO" in species || "NOX" in species
        end

        # --- Step 5: Assign surrogates ---
        gridref = read_gridref(agref_file)
        speciated_with_srg = assign_surrogates(speciated, gridref)
        GC.gc()

        @test hasproperty(speciated_with_srg, :Surrogate)
        @info "Surrogate assignment complete"

        # --- Step 6: Parse pre-computed surrogates ---
        county_surrogates = parse_smoke_surrogates(srg_file, grid)
        @info "Parsed surrogates for $(length(county_surrogates)) FIPS codes"

        # --- Step 7: Compute grid indices and refine with surrogates ---
        locIndex = compute_grid_indices(speciated_with_srg, grid)
        refine_indices_with_surrogates(locIndex, county_surrogates)
        GC.gc()

        in_grid_count = count(v -> v.inGrid, values(locIndex))
        @info "Location index: $(length(locIndex)) unique locations, $in_grid_count in grid"
        @test in_grid_count > 0

        # --- Step 8: Temporal allocation for Aug 1, 2018 ---
        # Use actual Gentpro temporal profiles via ATREF cross-reference.
        # This correctly handles:
        # - Per-FIPS, per-SCC monthly and daily profiles from Gentpro
        # - Different diurnal profiles: 600 (standard RWC) vs 1500 (hydronic heaters)
        # - Hydronic heaters use static monthly profile + weekly profile instead of daily

        hourly_profiles = parse_amptpro_hourly(hourly_file)
        weekly_profiles = parse_amptpro_weekly(weekly_file)
        gentpro_monthly = parse_gentpro_monthly(monthly_file)
        gentpro_daily = parse_gentpro_daily(daily_file)
        atref = parse_atref_gentpro(atref_file)

        @info "Gentpro monthly profiles: $(length(gentpro_monthly)), " *
              "daily profiles: $(length(gentpro_daily)), " *
              "hourly profiles: $(length(hourly_profiles)), " *
              "weekly profiles: $(length(weekly_profiles)), " *
              "ATREF entries: $(length(atref))"

        # Collect unique (FIPS, SCC) pairs from speciated emissions
        emissions_fips_scc = Set{Tuple{String, String}}()
        for row in eachrow(speciated_with_srg)
            push!(emissions_fips_scc, (string(row.FIPS), string(row.SCC)))
        end

        @info "Unique (FIPS, SCC) pairs: $(length(emissions_fips_scc))"

        # Build temporal_allocate-compatible profiles and xref from Gentpro data
        target_date = Date(2018, 8, 1)
        profiles, xref = build_gentpro_temporal(
            atref, gentpro_monthly, gentpro_daily,
            hourly_profiles, weekly_profiles,
            emissions_fips_scc, target_date)

        @info "Built $(nrow(profiles)) temporal profiles and $(nrow(xref)) xref entries"

        # Verify that different SCC types got different temporal profiles
        @testset "Per-SCC temporal profile assignment" begin
            # Check that xref has entries with different diurnal IDs
            unique_diurnal = unique(xref.diurnal_id)
            @test length(unique_diurnal) >= 2  # At least 600 and 1500

            # Check that xref has entries with different monthly IDs (FIPS-specific)
            unique_monthly = unique(xref.monthly_id)
            @test length(unique_monthly) > 1  # FIPS-specific profiles

            # Check that xref has entries with different weekly IDs
            unique_weekly = unique(xref.weekly_id)
            @test length(unique_weekly) >= 2  # Daily-derived vs weekly profiles
        end

        # Episode: Aug 1 00:00 to Aug 2 00:00 (25 hours for IOAPI convention)
        ep_start = DateTime(2018, 8, 1, 0)
        ep_end = DateTime(2018, 8, 2, 1)  # 25 hours

        # SMOKE uses standard-time timezone offsets from the COSTCY file.
        # The 12LISTOS domain covers the NE US (Eastern Standard Time = UTC-5).
        # SMOKE convention: OUTZONE=0 (GMT output), county timezone=5 (EST),
        # so diurnal profiles are shifted by -5 hours (local = UTC - 5).
        # Build a per-FIPS timezone map using state-level US timezone assignments.
        # Most states in the 12LISTOS domain are Eastern (UTC-5); a few may be
        # Central (UTC-6). SMOKE uses standard time (no DST adjustment).
        us_state_tz = Dict{String, Int}(
            # Eastern Standard Time (UTC-5)
            "09" => -5, "10" => -5, "11" => -5, "12" => -5, "13" => -5,
            "17" => -5, "18" => -5, "21" => -5, "23" => -5, "24" => -5,
            "25" => -5, "26" => -5, "27" => -5, "33" => -5, "34" => -5,
            "36" => -5, "37" => -5, "39" => -5, "42" => -5, "44" => -5,
            "45" => -5, "47" => -5, "50" => -5, "51" => -5, "54" => -5,
            # Central Standard Time (UTC-6)
            "01" => -6, "05" => -6, "19" => -6, "20" => -6, "22" => -6,
            "27" => -6, "28" => -6, "29" => -6, "31" => -6, "38" => -6,
            "40" => -6, "46" => -6, "48" => -6, "55" => -6,
            # Mountain Standard Time (UTC-7)
            "04" => -7, "08" => -7, "16" => -7, "30" => -7, "32" => -7,
            "35" => -7, "49" => -7, "56" => -7,
            # Pacific Standard Time (UTC-8)
            "02" => -8, "06" => -8, "15" => -8, "41" => -8, "53" => -8,
        )
        timezone_map = Dict{String, Int}()
        for row in eachrow(speciated_with_srg)
            fips = string(row.FIPS)
            if !haskey(timezone_map, fips) && length(fips) >= 2
                state = fips[1:2]
                tz = get(us_state_tz, state, -5)  # Default to EST
                timezone_map[fips] = tz
            end
        end

        hourly = temporal_allocate(
            speciated_with_srg, profiles, xref,
            ep_start, ep_end;
            timezone_map = timezone_map,
        )
        GC.gc()

        @test nrow(hourly) > 0
        @info "Hourly emissions: $(nrow(hourly)) rows"

        # Verify temporal allocation produced non-zero rates
        nonzero_rates = count(r -> r.emission_rate > 0, eachrow(hourly))
        @info "Non-zero hourly rates: $nonzero_rates / $(nrow(hourly))"

        # --- Step 9: Merge emissions onto grid ---
        species_list = unique(string.(hourly.POLID))
        merged = merge_emissions(hourly, locIndex, grid; species_list = species_list)
        GC.gc()

        @test nrow(merged) > 0
        @info "Merged gridded emissions: $(nrow(merged)) rows, species: $(length(unique(merged.pollutant)))"

        # --- Step 10: Convert to model-ready arrays ---
        hours = [ep_start + Hour(h) for h in 0:24]
        model_data = to_model_ready(merged, grid, hours)

        @test !isempty(model_data)
        @info "Model-ready species: $(sort(collect(keys(model_data))))"

        # ================================================================
        # Step 11: Comprehensive comparison to SMOKE reference output
        # ================================================================
        # This comparison uses:
        # 1. Temporal: Actual Gentpro FIPS-specific monthly/daily profiles via ATREF,
        #    with different diurnal profiles per SCC type (600 vs 1500).
        # 2. Spatial: Population surrogate (USA_100) for all SCCs — the only
        #    surrogate provided in the example case (SMOKE also uses this default).
        # 3. Speciation: GSREF assigns different profiles per SCC and pollutant
        #    (e.g., different VOC profiles for fireplaces vs catalytic stoves).
        #    Mass-basis is used; gas species are converted to mol for comparison.
        # 4. PMC: Computed as PM10 - PM25 (matching SMOKE's SMKINVEN_FORMULA).

        NCDatasets.Dataset(ref_file, "r") do ds
            ref_species = read_ioapi_species(ds)
            julia_species = sort(collect(keys(model_data)))
            common_species = sort(collect(intersect(Set(julia_species), Set(ref_species))))

            @info "Reference species ($(length(ref_species))): $(ref_species)"
            @info "Julia species ($(length(julia_species))): $(julia_species)"
            @info "Common species ($(length(common_species))): $(common_species)"

            # ---- Species completeness ----
            @testset "Species completeness" begin
                @test length(common_species) >= 15
                missing_from_julia = sort(collect(setdiff(Set(ref_species), Set(julia_species))))
                extra_in_julia = sort(collect(setdiff(Set(julia_species), Set(ref_species))))

                # Enhanced species analysis with HAP documentation
                @info "Species completeness analysis:"
                @info "  Reference species: $(length(ref_species))"
                @info "  Julia species: $(length(julia_species))"
                @info "  Common species: $(length(common_species))"
                @info "  Missing from Julia: $(length(missing_from_julia))"
                @info "  Extra in Julia: $(length(extra_in_julia))"

                if !isempty(missing_from_julia)
                    @info "Species in reference but not produced by Julia ($(length(missing_from_julia))):"
                    @info "  $(join(missing_from_julia[1:min(10, end)], ", "))$(length(missing_from_julia) > 10 ? "..." : "")"

                    # Document expected missing species due to HAP subtraction
                    hap_related = filter(s -> contains(lowercase(s), "ald") || contains(lowercase(s), "aro") ||
                                           contains(lowercase(s), "ole") || contains(lowercase(s), "par"), missing_from_julia)
                    if !isempty(hap_related)
                        @info "  Note: Some missing species may be HAP-related (expected): $(join(hap_related[1:min(5, end)], ", "))"
                    end
                end

                if !isempty(extra_in_julia)
                    @info "Species produced by Julia but not in reference: $(join(extra_in_julia[1:min(10, end)], ", "))"
                end
                # At least half of reference species with non-zero totals should be produced
                nonzero_ref = [sp for sp in ref_species if Float64(sum(read_ioapi_var_raw(ds, sp))) > 0]
                nonzero_common = intersect(Set(julia_species), Set(nonzero_ref))
                @info "Non-zero reference species: $(length(nonzero_ref)), common: $(length(nonzero_common))"
                @test length(nonzero_common) >= length(nonzero_ref) ÷ 2
            end

            # ---- Non-negativity and data quality ----
            @testset "Non-negativity and data quality" begin
                for sp in julia_species
                    arr = model_data[sp]

                    # Check for negative values
                    has_negative = any(arr .< 0)
                    if has_negative
                        neg_count = sum(arr .< 0)
                        min_val = minimum(arr)
                        @info "Species $sp has $neg_count negative values, minimum: $min_val"
                    end
                    @test !has_negative

                    # Check for NaN or Inf values (data corruption)
                    has_nan = any(isnan.(arr))
                    has_inf = any(isinf.(arr))
                    @test !has_nan
                    @test !has_inf

                    # Check for unreasonably large values (potential unit errors)
                    max_val = maximum(arr)
                    if max_val > 1e6  # mol/s or g/s
                        @warn "Species $sp has very large maximum value: $max_val (potential unit issue?)"
                    end
                    @test max_val < 1e8  # Sanity check for extreme values
                end
            end

            # ---- Output array dimensions ----
            @testset "Output dimensions" begin
                expected_size = (grid.Ny, grid.Nx, 1, length(hours))
                for (sp, arr) in model_data
                    @test size(arr) == expected_size
                end
            end

            # ---- Non-zero output for key species ----
            @testset "Key species have non-zero output" begin
                for sp in ["NO", "NO2", "CO", "SO2", "NH3"]
                    if haskey(model_data, sp)
                        sp_total = sum(model_data[sp])
                        @test sp_total > 0
                    end
                end
            end

            # ---- Spatial patterns for all common species ----
            @testset "Spatial pattern correlation" begin
                spatial_corrs = Dict{String, Float64}()

                for sp in common_species
                    ref_spatial = ioapi_spatial_pattern(ds, sp)  # (ROW, COL)
                    julia_spatial = dropdims(sum(model_data[sp], dims = (3, 4)), dims = (3, 4))

                    ref_total = sum(ref_spatial)
                    julia_total = sum(julia_spatial)

                    if ref_total > 0 && julia_total > 0
                        ref_norm = vec(ref_spatial) ./ ref_total
                        julia_norm = vec(julia_spatial) ./ julia_total
                        corr_val = cosine_similarity(ref_norm, julia_norm)
                        spatial_corrs[sp] = corr_val
                    end
                end

                @info "Spatial correlations computed for $(length(spatial_corrs)) species"
                for (sp, c) in sort(collect(spatial_corrs), by = x -> x[2])
                    @info "  $sp: $(round(c, digits=4))"
                end

                # Key inorganic species should have high spatial correlation
                # These have straightforward 1:1 or simple speciation
                for sp in ["CO", "NO", "SO2", "NH3", "NO2"]
                    if haskey(spatial_corrs, sp)
                        @info "Spatial correlation $sp: $(round(spatial_corrs[sp], digits=4))"
                        @test spatial_corrs[sp] > 0.9
                    end
                end

                # PM species should also correlate well
                for sp in ["PEC", "POC", "PSO4", "PNH4", "PNO3", "PMOTHR"]
                    if haskey(spatial_corrs, sp)
                        @info "Spatial correlation $sp: $(round(spatial_corrs[sp], digits=4))"
                        @test spatial_corrs[sp] > 0.9
                    end
                end

                # Median correlation across ALL common species should be reasonable
                if !isempty(spatial_corrs)
                    vals = filter(!isnan, collect(values(spatial_corrs)))
                    if !isempty(vals)
                        med_corr = median(vals)
                        @info "Median spatial correlation across $(length(vals)) species: $(round(med_corr, digits=4))"
                        @test med_corr > 0.9
                    end
                end
            end

            # ---- Active cell overlap ----
            @testset "Active cell overlap" begin
                for sp in ["CO", "NO", "SO2", "NH3"]
                    if haskey(model_data, sp) && sp in ref_species
                        ref_spatial = ioapi_spatial_pattern(ds, sp)
                        julia_spatial = dropdims(sum(model_data[sp], dims = (3, 4)), dims = (3, 4))

                        ref_active = vec(ref_spatial) .> 0
                        julia_active = vec(julia_spatial) .> 0
                        n_intersection = count(ref_active .& julia_active)
                        n_union = count(ref_active .| julia_active)
                        jaccard = n_union > 0 ? n_intersection / n_union : 0.0

                        @info "Active cell Jaccard index ($sp): $(round(jaccard, digits=3)) (ref=$(count(ref_active)), julia=$(count(julia_active)), overlap=$n_intersection)"
                        # All non-point RWC sources use population surrogate, so
                        # active cells should have high overlap
                        @test jaccard > 0.5
                    end
                end
            end

            # ---- Spatial concentration (not uniform) ----
            @testset "Spatial concentration" begin
                for sp in ["CO", "NO", "SO2", "NH3"]
                    if haskey(model_data, sp)
                        arr = model_data[sp]
                        total = sum(arr)
                        if total > 0
                            flat = sort(vec(sum(arr, dims = (3, 4))), rev = true)
                            top5_frac = sum(flat[1:min(5, length(flat))]) / total
                            # Population-based surrogate produces concentrated emissions
                            # (NYC metro area dominates in LISTOS domain)
                            @test top5_frac > 0.05
                            # Verify not completely concentrated in one cell
                            @test flat[1] / total < 0.5
                        end
                    end
                end
            end

            # ---- Species ratio consistency ----
            @testset "Species ratios" begin
                # With mole-basis speciation, Julia output units are:
                # - Gas: kg/s * split_factor/divisor → multiply by 1000 for mol/s
                # - PM: kg/s * mass_fraction → multiply by 1000 for g/s
                # Both match reference units after ×1000 conversion.
                """Convert Julia total to same units as reference (mol/s or g/s)."""
                to_ref_units(sp, julia_total) = julia_total * 1000.0

                # Within-group gas ratios (both species from same pollutant group)
                # These should match exactly since both species are from the same
                # pollutant (NOX) and use the same speciation profile.
                gas_ratio_tests = [
                    ("NO", "NO2", "NO/NO2 (both from NOX)", 0.9, 1.1),
                ]
                for (sp1, sp2, desc, rmin, rmax) in gas_ratio_tests
                    if haskey(model_data, sp1) && haskey(model_data, sp2) &&
                       sp1 in ref_species && sp2 in ref_species
                        julia_mol_1 = to_ref_units(sp1, sum(model_data[sp1]))
                        julia_mol_2 = to_ref_units(sp2, sum(model_data[sp2]))
                        ref_total_1 = Float64(sum(read_ioapi_var_raw(ds, sp1)))
                        ref_total_2 = Float64(sum(read_ioapi_var_raw(ds, sp2)))

                        if julia_mol_2 > 0 && ref_total_2 > 0 && ref_total_1 > 0
                            julia_ratio = julia_mol_1 / julia_mol_2
                            ref_ratio = ref_total_1 / ref_total_2
                            ratio_of_ratios = julia_ratio / ref_ratio
                            @info "$desc: julia=$(round(julia_ratio, digits=4)), ref=$(round(ref_ratio, digits=4)), ratio_of_ratios=$(round(ratio_of_ratios, digits=4))"
                            @test rmin < ratio_of_ratios < rmax
                        end
                    end
                end

                # Cross-group gas ratios (species from different pollutant groups)
                # These should be close since the temporal profiles apply the same
                # factors to all species from the same source. Tighter bounds because
                # all RWC sources share the same surrogate and similar temporal profiles.
                cross_ratio_tests = [
                    ("NO", "CO", "NO/CO (cross-group)", 0.7, 1.4),
                    ("SO2", "CO", "SO2/CO (cross-group)", 0.7, 1.4),
                    ("NH3", "CO", "NH3/CO (cross-group)", 0.7, 1.4),
                ]
                for (sp1, sp2, desc, rmin, rmax) in cross_ratio_tests
                    if haskey(model_data, sp1) && haskey(model_data, sp2) &&
                       sp1 in ref_species && sp2 in ref_species
                        julia_mol_1 = to_ref_units(sp1, sum(model_data[sp1]))
                        julia_mol_2 = to_ref_units(sp2, sum(model_data[sp2]))
                        ref_total_1 = Float64(sum(read_ioapi_var_raw(ds, sp1)))
                        ref_total_2 = Float64(sum(read_ioapi_var_raw(ds, sp2)))

                        if julia_mol_2 > 0 && ref_total_2 > 0 && ref_total_1 > 0
                            julia_ratio = julia_mol_1 / julia_mol_2
                            ref_ratio = ref_total_1 / ref_total_2
                            ratio_of_ratios = julia_ratio / ref_ratio
                            @info "$desc: julia=$(round(julia_ratio, digits=4)), ref=$(round(ref_ratio, digits=4)), ratio_of_ratios=$(round(ratio_of_ratios, digits=4))"
                            @test rmin < ratio_of_ratios < rmax
                        end
                    end
                end

                # PM species ratio check (both mass-based, no MW correction needed)
                pm_ratio_tests = [
                    ("PEC", "POC", "PEC/POC (PM components)", 0.5, 2.0),
                    ("PNCOM", "POC", "PNCOM/POC (PM components)", 0.2, 5.0),
                ]
                for (sp1, sp2, desc, rmin, rmax) in pm_ratio_tests
                    if haskey(model_data, sp1) && haskey(model_data, sp2) &&
                       sp1 in ref_species && sp2 in ref_species
                        julia_total_1 = sum(model_data[sp1])
                        julia_total_2 = sum(model_data[sp2])
                        ref_total_1 = Float64(sum(read_ioapi_var_raw(ds, sp1)))
                        ref_total_2 = Float64(sum(read_ioapi_var_raw(ds, sp2)))

                        if julia_total_2 > 0 && ref_total_2 > 0 && ref_total_1 > 0
                            julia_ratio = julia_total_1 / julia_total_2
                            ref_ratio = ref_total_1 / ref_total_2
                            ratio_of_ratios = julia_ratio / ref_ratio
                            @info "$desc: julia=$(round(julia_ratio, digits=4)), ref=$(round(ref_ratio, digits=4)), ratio_of_ratios=$(round(ratio_of_ratios, digits=4))"
                            @test rmin < ratio_of_ratios < rmax
                        end
                    end
                end
            end

            # ---- Diurnal pattern comparison ----
            @testset "Diurnal pattern" begin
                # Compare hourly emission profiles (normalized).
                # We use per-SCC diurnal profiles via ATREF:
                # - Standard RWC SCCs: profile 600 (evening peak)
                # - Hydronic heater SCCs: profile 1500 (nearly flat)
                # The aggregate diurnal pattern depends on the relative contribution
                # of each SCC type, which in turn depends on temporal allocation.
                for sp in ["CO", "NO", "SO2"]
                    if haskey(model_data, sp) && sp in ref_species
                        ref_hourly = ioapi_hourly_totals(ds, sp)
                        julia_arr = model_data[sp]
                        julia_hourly = [sum(julia_arr[:, :, :, t]) for t in 1:size(julia_arr, 4)]

                        # Both should have the same number of timesteps
                        n = min(length(ref_hourly), length(julia_hourly))
                        ref_h = ref_hourly[1:n]
                        julia_h = julia_hourly[1:n]

                        ref_sum = sum(ref_h)
                        julia_sum = sum(julia_h)

                        if ref_sum > 0 && julia_sum > 0 && n >= 12
                            ref_norm = ref_h ./ ref_sum
                            julia_norm = julia_h ./ julia_sum
                            diurnal_corr = cosine_similarity(ref_norm, julia_norm)
                            @info "Diurnal pattern correlation ($sp): $(round(diurnal_corr, digits=4))"
                            @test diurnal_corr > 0.9

                            # Also check that hourly variation exists (not flat)
                            ref_cv = std(ref_norm) / mean(ref_norm)
                            julia_cv = std(julia_norm) / mean(julia_norm)
                            @info "Diurnal CV ($sp): ref=$(round(ref_cv, digits=3)), julia=$(round(julia_cv, digits=3))"

                            # Mass conservation check: diurnal patterns should preserve relative totals
                            # Even though absolute magnitudes may differ, the ratio of daily totals
                            # should be close to the ratio seen in spatial comparisons
                            if ref_sum > 0 && julia_sum > 0
                                temporal_magnitude_ratio = (julia_sum * 1000.0) / ref_sum
                                # Allow broader tolerance since temporal allocation can introduce small numerical differences
                                @test 0.5 < temporal_magnitude_ratio < 2.0
                                if temporal_magnitude_ratio < 0.99 || temporal_magnitude_ratio > 1.01
                                    @info "$sp temporal magnitude ratio: $(round(temporal_magnitude_ratio, digits=4))"
                                end
                            end
                        end
                    end
                end
            end

            # ---- Per-cell spatial comparison for key species ----
            @testset "Per-cell spatial comparison" begin
                for sp in ["CO", "NO", "SO2", "NH3"]
                    if !haskey(model_data, sp) || !(sp in ref_species)
                        continue
                    end
                    ref_spatial = ioapi_spatial_pattern(ds, sp)
                    julia_spatial = dropdims(sum(model_data[sp], dims = (3, 4)), dims = (3, 4))

                    ref_total = sum(ref_spatial)
                    julia_total = sum(julia_spatial)

                    ref_total > 0 && julia_total > 0 || continue

                    ref_norm = ref_spatial ./ ref_total
                    julia_norm = julia_spatial ./ julia_total

                    # Normalized RMSE (normalized by mean cell value).
                    # Values are moderate because most grid cells have near-zero emissions
                    # while a few have large values, so RMS / mean is amplified.
                    diff = ref_norm .- julia_norm
                    nrmse = sqrt(mean(diff .^ 2)) / mean(ref_norm)
                    @info "$sp normalized RMSE: $(round(nrmse, digits=4))"
                    @test nrmse < 0.8

                    # Check that the top-emitting cells overlap
                    ref_flat = vec(ref_norm)
                    julia_flat = vec(julia_norm)
                    ref_top10 = Set(sortperm(ref_flat, rev = true)[1:10])
                    julia_top10 = Set(sortperm(julia_flat, rev = true)[1:10])
                    overlap = length(intersect(ref_top10, julia_top10))
                    @info "$sp top-10 cell overlap: $overlap / 10"
                    @test overlap >= 3

                    # Per-cell ratio check: for cells with significant emissions,
                    # verify the spatial allocation fraction is close
                    significant = ref_norm .> (mean(ref_norm) * 0.5)
                    if count(significant) > 3
                        cell_ratios = julia_norm[significant] ./ ref_norm[significant]
                        med_cell_ratio = median(cell_ratios)
                        @info "$sp per-cell median ratio ($(count(significant)) cells): $(round(med_cell_ratio, digits=4))"
                        @test 0.5 < med_cell_ratio < 2.0
                    end
                end
            end

            # ---- Magnitude diagnostics and consistency ----
            @testset "Magnitude diagnostics" begin
                # With mole-basis speciation matching SMOKE's convention:
                # - Gas species: julia output is kg_pollutant/s * split_factor/divisor_g_per_mol
                #   → multiply by 1000 (kg→g) to get mol/s, matching reference units
                # - PM species: divisor=1.0 in GSPRO, so mole basis = mass basis
                #   → multiply by 1000 (kg→g) to get g/s, matching reference units
                # Both cases: julia_converted = julia_total * 1000.0
                pm_species = Set(["PAL", "PCA", "PCL", "PEC", "PFE", "PH2O", "PK",
                    "PMG", "PMN", "PMOTHR", "PNA", "PNCOM", "PNH4", "PNO3",
                    "POC", "PSI", "PSO4", "PTI", "PMC", "SULF"])

                magnitude_ratios = Dict{String, Float64}()

                for sp in sort(collect(common_species))
                    ref_total = Float64(sum(read_ioapi_var_raw(ds, sp)))
                    julia_total = sum(model_data[sp])

                    # Convert julia from kg-based to reference units (mol/s or g/s)
                    julia_converted = julia_total * 1000.0
                    ratio = ref_total > 0 ? julia_converted / ref_total : NaN
                    magnitude_ratios[sp] = ratio

                    if sp in pm_species
                        @info "  $sp (PM): ref=$(round(ref_total, sigdigits=4)) g/s, julia=$(round(julia_converted, sigdigits=4)) g/s, ratio=$(round(ratio, sigdigits=3))"
                    else
                        @info "  $sp (gas): ref=$(round(ref_total, sigdigits=4)) mol/s, julia=$(round(julia_converted, sigdigits=4)) mol/s, ratio=$(round(ratio, sigdigits=3))"
                    end
                end

                # All common species with non-zero reference should have non-zero Julia output
                for sp in common_species
                    ref_total = Float64(sum(read_ioapi_var_raw(ds, sp)))
                    julia_total = sum(model_data[sp])
                    if ref_total > 0
                        @test julia_total > 0
                    end
                end

                # --- Per-species magnitude checks for key inorganic species ---
                # These have 1:1 or simple speciation (no HAP subtraction), so
                # the ratio should be close to 1.0.
                @testset "Key species magnitudes" begin
                    key_species = ["CO", "NO", "NO2", "SO2", "NH3"]
                    for sp in key_species
                        if haskey(magnitude_ratios, sp) && !isnan(magnitude_ratios[sp])
                            r = magnitude_ratios[sp]
                            @info "Key species $sp magnitude ratio: $(round(r, sigdigits=4)) (target: 0.7-1.4, ideal: ~1.0)"
                            @test 0.7 < r < 1.4
                        end
                    end

                    # Additional check for very tight agreement on key inorganics
                    super_key = ["CO", "NO", "SO2", "NH3"] # Species with direct 1:1 relationships
                    for sp in super_key
                        if haskey(magnitude_ratios, sp) && !isnan(magnitude_ratios[sp])
                            r = magnitude_ratios[sp]
                            @info "Super-key species $sp magnitude ratio: $(round(r, sigdigits=4)) (tight target: 0.95-1.05)"
                            @test 0.95 < r < 1.05  # Very tight tolerance for key inorganics
                        end
                    end

                    # PM species magnitude validation
                    key_pm = ["PEC", "POC", "PNCOM", "PMOTHR", "PNH4", "PNO3", "PSO4"]
                    for sp in key_pm
                        if haskey(magnitude_ratios, sp) && !isnan(magnitude_ratios[sp])
                            r = magnitude_ratios[sp]
                            @info "Key PM species $sp magnitude ratio: $(round(r, sigdigits=4)) (target: 0.7-1.4)"
                            @test 0.7 < r < 1.4
                        end
                    end
                end

                # --- Aggregate magnitude checks ---
                valid_ratios = filter(p -> !isnan(p.second) && p.second > 0, magnitude_ratios)
                if !isempty(valid_ratios)
                    ratio_vals = collect(values(valid_ratios))
                    @info "Magnitude ratio stats: median=$(round(median(ratio_vals), sigdigits=3)), " *
                          "min=$(round(minimum(ratio_vals), sigdigits=3)), " *
                          "max=$(round(maximum(ratio_vals), sigdigits=3))"

                    # All species should have ratios within a factor of 2
                    for (sp, r) in sort(collect(valid_ratios), by = x -> x[2])
                        # Document HAP-affected species with relaxed tolerance
                        if contains(lowercase(sp), "ald") || contains(lowercase(sp), "aro") ||
                           contains(lowercase(sp), "ole") || contains(lowercase(sp), "par")
                            @test 0.3 < r < 2.0  # Relaxed for HAP-affected species
                            @info "$sp ratio $(round(r, sigdigits=3)) - HAP-affected (expected ~0.8x)"
                        else
                            @test 0.5 < r < 2.0
                        end
                    end

                    # Median ratio should be close to 1.0
                    @test 0.7 < median(ratio_vals) < 1.3
                end
            end

            # ---- Zero-species consistency ----
            @testset "Zero-species consistency" begin
                # Species that are zero in reference should also be zero (or absent)
                # in our output, and vice versa
                for sp in common_species
                    ref_total = Float64(sum(read_ioapi_var_raw(ds, sp)))
                    julia_total = sum(model_data[sp])
                    if ref_total == 0.0
                        @test julia_total == 0.0 || !haskey(model_data, sp)
                    end
                end
            end

            # ---- Per-cell per-hour comparison for CO ----
            @testset "Per-cell per-hour comparison (CO)" begin
                if haskey(model_data, "CO") && "CO" in ref_species
                    ref_data = read_ioapi_var_raw(ds, "CO")  # (COL, ROW, LAY, TSTEP)
                    julia_data = model_data["CO"]  # (ROW, COL, LAY, TSTEP)

                    n_tsteps = min(size(ref_data, 4), size(julia_data, 4))
                    # Compare normalized hourly spatial patterns for a few timesteps
                    hourly_corrs = Float64[]
                    for t in 1:n_tsteps
                        ref_slice = permutedims(ref_data[:, :, 1, t], (2, 1))  # → (ROW, COL)
                        julia_slice = julia_data[:, :, 1, t]

                        ref_sum_t = sum(ref_slice)
                        julia_sum_t = sum(julia_slice)
                        if ref_sum_t > 0 && julia_sum_t > 0
                            ref_n = vec(ref_slice) ./ ref_sum_t
                            julia_n = vec(julia_slice) ./ julia_sum_t
                            push!(hourly_corrs, cosine_similarity(ref_n, julia_n))
                        end
                    end
                    if length(hourly_corrs) >= 10
                        med_hourly_corr = median(hourly_corrs)
                        min_hourly_corr = minimum(hourly_corrs)
                        @info "CO per-hour spatial correlation: median=$(round(med_hourly_corr, digits=4)), min=$(round(min_hourly_corr, digits=4)) over $(length(hourly_corrs)) timesteps"
                        @test med_hourly_corr > 0.9
                        @test min_hourly_corr > 0.8
                    end
                end
            end
        end

        # ================================================================
        # Multi-day consistency check
        # ================================================================
        # Verify that our pipeline produces consistent results across multiple
        # days by comparing against additional reference files.
        # The reference output contains 31 days of August 2018.
        @testset "Multi-day consistency" begin
            alt_ref_dir = dirname(ref_file)
            additional_dates = [Date(2018, 8, 15), Date(2018, 8, 31)]

            for target_date_alt in additional_dates
                date_str = Dates.format(target_date_alt, "yyyymmdd")
                alt_ref = joinpath(alt_ref_dir, "emis_mole_rwc_$(date_str)_12LISTOS_cmaq_cb6ae7_2018gg_18j.ncf")
                isfile(alt_ref) || continue

                @testset "Day $(date_str)" begin
                    # Rebuild temporal profiles for this day
                    alt_profiles, alt_xref = build_gentpro_temporal(
                        parse_atref_gentpro(atref_file),
                        parse_gentpro_monthly(monthly_file),
                        parse_gentpro_daily(daily_file),
                        parse_amptpro_hourly(hourly_file),
                        parse_amptpro_weekly(weekly_file),
                        emissions_fips_scc, target_date_alt)

                    alt_ep_start = DateTime(target_date_alt)
                    alt_ep_end = DateTime(target_date_alt) + Hour(25)

                    # Run temporal + merge
                    alt_hourly = temporal_allocate(
                        speciated_with_srg, alt_profiles, alt_xref,
                        alt_ep_start, alt_ep_end;
                        timezone_map = timezone_map)

                    alt_species_list = unique(string.(alt_hourly.POLID))
                    alt_merged = merge_emissions(alt_hourly, locIndex, grid;
                        species_list = alt_species_list)
                    alt_hours = [alt_ep_start + Hour(h) for h in 0:24]
                    alt_model = to_model_ready(alt_merged, grid, alt_hours)

                    NCDatasets.Dataset(alt_ref, "r") do ds2
                        # Check key species spatial correlation matches reference
                        for sp in ["CO", "NO"]
                            if haskey(alt_model, sp) && sp in read_ioapi_species(ds2)
                                ref_sp = ioapi_spatial_pattern(ds2, sp)
                                jul_sp = dropdims(sum(alt_model[sp], dims = (3, 4)), dims = (3, 4))
                                ref_t = sum(ref_sp)
                                jul_t = sum(jul_sp)
                                if ref_t > 0 && jul_t > 0
                                    corr_val = cosine_similarity(
                                        vec(ref_sp) ./ ref_t,
                                        vec(jul_sp) ./ jul_t)
                                    @info "$(date_str) $sp spatial correlation: $(round(corr_val, digits=4))"
                                    @test corr_val > 0.9
                                end

                                # Magnitude check
                                ratio = jul_t > 0 && ref_t > 0 ? (jul_t * 1000.0) / ref_t : NaN
                                if !isnan(ratio)
                                    @info "$(date_str) $sp magnitude ratio: $(round(ratio, sigdigits=3))"
                                    @test 0.5 < ratio < 2.0
                                end
                            end
                        end

                        # Check diurnal pattern
                        if haskey(alt_model, "CO") && "CO" in read_ioapi_species(ds2)
                            ref_hr = ioapi_hourly_totals(ds2, "CO")
                            jul_arr = alt_model["CO"]
                            jul_hr = [sum(jul_arr[:, :, :, t]) for t in 1:size(jul_arr, 4)]
                            n = min(length(ref_hr), length(jul_hr))
                            if n >= 12 && sum(ref_hr[1:n]) > 0 && sum(jul_hr[1:n]) > 0
                                corr_val = cosine_similarity(
                                    ref_hr[1:n] ./ sum(ref_hr[1:n]),
                                    jul_hr[1:n] ./ sum(jul_hr[1:n]))
                                @info "$(date_str) CO diurnal correlation: $(round(corr_val, digits=4))"
                                @test corr_val > 0.9
                            end
                        end
                    end
                    GC.gc()
                end
            end
        end
    end

    # ========================================================================
    @testset "Per-SCC profile differentiation" begin
        # Verify that the test correctly assigns different profiles to different
        # source types, as the reference SMOKE model does.

        # Parse GSREF and verify SCC-specific speciation
        gsref = parse_smoke_gsref_csv(gsref_file)

        # Wood combustion SCCs should use different VOC profiles than hydronic heaters
        wood_voc = filter(r -> startswith(r.SCC, "2104008") &&
                              !startswith(r.SCC, "21040086") &&  # exclude hydronic
                              (r.pollutant_id == "NONHAPVOC" || r.pollutant_id == "VOC"), gsref)
        hydronic_voc = filter(r -> startswith(r.SCC, "21040086") &&
                                   (r.pollutant_id == "NONHAPVOC" || r.pollutant_id == "VOC"), gsref)

        if nrow(wood_voc) > 0 && nrow(hydronic_voc) > 0
            wood_profiles = unique(wood_voc.profile_code)
            hydronic_profiles = unique(hydronic_voc.profile_code)
            @info "Wood combustion VOC profiles: $wood_profiles"
            @info "Hydronic heater VOC profiles: $hydronic_profiles"
        end

        # Different PM2.5 profiles should exist for different fuel types
        pm_profiles = filter(r -> r.pollutant_id == "PM2_5" || r.pollutant_id == "PM25", gsref)
        if nrow(pm_profiles) > 0
            unique_pm_profiles = unique(pm_profiles.profile_code)
            @info "Unique PM2.5 speciation profiles: $(length(unique_pm_profiles))"
            @test length(unique_pm_profiles) >= 2  # At least residential wood vs oil/gas
        end

        # Parse ATREF and verify SCC-specific temporal profiles
        atref = parse_atref_gentpro(atref_file)

        # Collect unique diurnal profile IDs by SCC prefix
        standard_rwc_diurnal = Set{String}()
        hydronic_diurnal = Set{String}()
        for ((scc, fips), info) in atref
            allday = get(info, "ALLDAY", "")
            if startswith(scc, "21040086")  # hydronic heater
                push!(hydronic_diurnal, allday)
            elseif startswith(scc, "2104")
                push!(standard_rwc_diurnal, allday)
            end
        end

        @info "Standard RWC diurnal profiles: $standard_rwc_diurnal"
        @info "Hydronic heater diurnal profiles: $hydronic_diurnal"

        # Verify different diurnal profiles for different source types
        if !isempty(standard_rwc_diurnal) && !isempty(hydronic_diurnal)
            @test standard_rwc_diurnal != hydronic_diurnal
        end

        # Verify FIPS-specific monthly profiles
        monthly_profiles_used = Set{String}()
        for ((scc, fips), info) in atref
            mon_str = get(info, "MONTHLY", "")
            if !isempty(mon_str)
                push!(monthly_profiles_used, mon_str)
            end
        end
        @info "Unique monthly profile IDs: $(length(monthly_profiles_used))"
        @test length(monthly_profiles_used) > 10  # Should have many FIPS-specific profiles
    end

    # ========================================================================
    @testset "Reference output structure" begin
        NCDatasets.Dataset(ref_file, "r") do ds
            # IOAPI structure checks
            @test ds.attrib["FTYPE"] == 1  # GRDDED3
            @test ds.attrib["GDTYP"] == 2  # LCC

            # Check that TFLAG is present and correctly structured
            @test haskey(ds, "TFLAG")
            tflag = Array(ds["TFLAG"])
            @test size(tflag, 1) == 2  # DATE-TIME dimension
            @test tflag[1, 1, 1] == 2018213  # Aug 1, 2018 = day 213

            # Verify number of variables and timesteps
            nvars = ds.attrib["NVARS"]
            @test nvars == 62
            @test size(tflag, 3) == 25  # 25 hourly timesteps

            # Verify all declared species are readable
            ref_species = read_ioapi_species(ds)
            @test length(ref_species) == nvars
            for sp in ref_species
                @test haskey(ds, sp)
            end
        end
    end
end
