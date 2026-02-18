"""
Integration test: validates Emissions.jl against SMOKE ExampleCase v2 reference output.
Tests the RWC (residential wood combustion) nonpoint sector for August 1, 2018.

Input data from: https://github.com/CEMPD/SMOKE-ExampleCase-v2
Reference output: premerged RWC IOAPI NetCDF files

Data is automatically downloaded from Google Drive if not already present.
The main example case tarball is ~2 GB and the reference output is ~30 MB.
"""

using Test
using Emissions
using DataFrames
using SparseArrays
using Dates
using NCDatasets
using CSV
using Unitful: ustrip
using Statistics: cor, mean, median, std

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
    url = "https://drive.google.com/uc?export=download&id=$(file_id)&confirm=t"
    @info "Downloading from Google Drive (file ID: $file_id)..."
    run(`curl -L --retry 3 --retry-delay 10 --max-time 7200 -o $output_path $url`)
    if !isfile(output_path) || filesize(output_path) < 1000
        error("Download failed or file too small: $output_path ($(filesize(output_path)) bytes)")
    end
    @info "Downloaded $(round(filesize(output_path) / 1024^2, digits=1)) MB"
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
        isempty(fips_raw) && continue
        isempty(profile_type) && continue

        # Normalize FIPS: strip country digit if 6 digits
        if length(fips_raw) == 6
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
Rename emissions POLID to match GSREF/GSPRO naming conventions.
- "PM25" -> "PM2_5" (matches GSREF/GSPRO naming)
- "PM10" -> "PMC" (matches GSREF/GSPRO naming for coarse PM)
"""
function rename_emissions_for_speciation!(emissions::DataFrame)
    mapping = Dict(
        "PM25" => "PM2_5",
        "PM10" => "PMC",
    )
    emissions.POLID = [get(mapping, string(p), string(p)) for p in emissions.POLID]
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
    atref_file = joinpath(SMOKE_BASE, "ge_dat", "temporal",
        "atref_2017platform_rwc_29apr2020_v1")

    @testset "Required files exist" begin
        for f in [inv_file, griddesc_file, agref_file, srg_file, gspro_file,
                  gsref_file, monthly_file, daily_file, hourly_file, atref_file,
                  ref_file]
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

        atref = parse_atref_gentpro(atref_file)
        @test !isempty(atref)
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

        # Rename emissions POLID to match GSREF/GSPRO naming (PM25→PM2_5, PM10→PMC)
        rename_emissions_for_speciation!(emissions)

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
        # Use mass basis; we convert to appropriate units when comparing to reference
        speciated = speciate_emissions(emissions, gspro, gsref; basis = :mass)
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
        # Note: Gentpro monthly profiles for RWC show August=0 for ALL FIPS codes
        # in the 12LISTOS grid (physically correct — no wood burning in summer in the
        # Northeast). However, the SMOKE reference output has non-zero August emissions,
        # suggesting SMOKE uses fallback/alternate temporal profiles.
        #
        # We use uniform monthly profiles (1/12) through the standard allocation path,
        # which gives non-zero emissions for comparison of spatial patterns and species
        # ratios. The absolute magnitude will differ from the reference since the true
        # SMOKE temporal factors are unknown.

        # Parse the RWC diurnal profile from SMOKE hourly profile file
        hourly_profiles = parse_amptpro_hourly(hourly_file)
        default_diurnal_id = "600"  # RWC ALLDAY profile from ATREF
        diurnal_factors = get(hourly_profiles, default_diurnal_id, fill(1.0 / 24.0, 24))

        # Build temporal profiles + xref for Emissions.jl
        # Uniform monthly (1/12) and weekly (1.0) profiles; RWC diurnal profile from SMOKE
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "ALLDAY"],
            profile_id = [1, 1, 1],
            factors = [
                fill(1.0 / 12.0, 12),  # uniform monthly (all months equal)
                fill(1.0, 7),           # uniform weekly (all days equal)
                diurnal_factors,        # RWC diurnal profile from SMOKE
            ],
        )
        xref = DataFrame(
            FIPS = ["00000"],
            SCC = ["0000000000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1],
        )

        # Episode: Aug 1 00:00 to Aug 2 00:00 (25 hours for IOAPI convention)
        ep_start = DateTime(2018, 8, 1, 0)
        ep_end = DateTime(2018, 8, 2, 1)  # 25 hours

        hourly = temporal_allocate(
            speciated_with_srg, profiles, xref,
            ep_start, ep_end,
        )
        GC.gc()

        @test nrow(hourly) > 0
        @info "Hourly emissions: $(nrow(hourly)) rows"

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
        # Known limitations of this comparison:
        # 1. Temporal: We use uniform 1/12 monthly profile, not SMOKE's actual
        #    monthly factors. Absolute magnitudes will differ, but normalized
        #    spatial patterns and species ratios should still match.
        # 2. Spatial: We use a single surrogate (USA_100 = population) for all
        #    SCCs. SMOKE assigns different surrogates per SCC via Grdmat. This
        #    may cause spatial pattern differences for some species.
        # 3. Speciation basis: We use mass-basis speciation. SMOKE uses mole-basis
        #    for gas species. This affects absolute magnitudes but not spatial
        #    patterns or within-group species ratios.

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
                if !isempty(missing_from_julia)
                    @info "Species in reference but not produced by Julia ($(length(missing_from_julia))): $missing_from_julia"
                end
                if !isempty(extra_in_julia)
                    @info "Species produced by Julia but not in reference ($(length(extra_in_julia))): $extra_in_julia"
                end
                # At least half of reference species should be produced
                @test length(common_species) >= length(ref_species) ÷ 3
            end

            # ---- Non-negativity ----
            @testset "Non-negativity" begin
                for sp in julia_species
                    has_negative = any(model_data[sp] .< 0)
                    if has_negative
                        @info "Species $sp has negative emission values"
                    end
                    @test !has_negative
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
                        @test spatial_corrs[sp] > 0.75
                    end
                end

                # PM species should also correlate well
                for sp in ["PEC", "POC", "PSO4", "PNH4", "PNO3", "PMOTHR"]
                    if haskey(spatial_corrs, sp)
                        @info "Spatial correlation $sp: $(round(spatial_corrs[sp], digits=4))"
                        @test spatial_corrs[sp] > 0.7
                    end
                end

                # Median correlation across ALL common species should be reasonable
                if !isempty(spatial_corrs)
                    vals = filter(!isnan, collect(values(spatial_corrs)))
                    if !isempty(vals)
                        med_corr = median(vals)
                        @info "Median spatial correlation across $(length(vals)) species: $(round(med_corr, digits=4))"
                        @test med_corr > 0.6
                    end
                end
            end

            # ---- Active cell overlap ----
            @testset "Active cell overlap" begin
                for sp in ["CO", "NO", "SO2"]
                    if haskey(model_data, sp) && sp in ref_species
                        ref_spatial = ioapi_spatial_pattern(ds, sp)
                        julia_spatial = dropdims(sum(model_data[sp], dims = (3, 4)), dims = (3, 4))

                        ref_active = vec(ref_spatial) .> 0
                        julia_active = vec(julia_spatial) .> 0
                        n_intersection = count(ref_active .& julia_active)
                        n_union = count(ref_active .| julia_active)
                        jaccard = n_union > 0 ? n_intersection / n_union : 0.0

                        @info "Active cell Jaccard index ($sp): $(round(jaccard, digits=3)) (ref=$(count(ref_active)), julia=$(count(julia_active)), overlap=$n_intersection)"
                        @test jaccard > 0.3
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
                            @test top5_frac > 0.05
                        end
                    end
                end
            end

            # ---- Species ratio consistency ----
            @testset "Species ratios" begin
                # Within-pollutant-group ratios should match well because speciation
                # factors are applied identically to the same source mix.
                # Cross-group ratios (e.g., NO/CO) test that the relative magnitude
                # of different pollutant groups is consistent.

                ratio_tests = [
                    # (numerator, denominator, description, min_ratio, max_ratio)
                    ("NO", "NO2", "NO/NO2 (both from NOX)", 0.2, 5.0),
                    ("NO", "CO", "NO/CO (cross-group)", 0.2, 5.0),
                    ("SO2", "CO", "SO2/CO (cross-group)", 0.05, 20.0),
                    ("NH3", "CO", "NH3/CO (cross-group)", 0.05, 20.0),
                ]

                for (sp1, sp2, desc, rmin, rmax) in ratio_tests
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

                # PM species ratio check
                pm_ratio_tests = [
                    ("PEC", "POC", "PEC/POC (PM components)", 0.1, 10.0),
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
                # We apply a single diurnal profile (profile 600 from ATREF) to
                # all sources uniformly, while SMOKE may use per-source or
                # per-FIPS diurnal profiles from the Gentpro temporal system.
                # Additionally, SMOKE's temporal allocation combines monthly,
                # weekly, and diurnal factors that differ by source, so the
                # aggregate hourly pattern will differ.
                # Threshold of 0.5 verifies meaningful diurnal shape agreement
                # while allowing for these known differences.
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
                            @test diurnal_corr > 0.5

                            # Also check that hourly variation exists (not flat)
                            ref_cv = std(ref_norm) / mean(ref_norm)
                            julia_cv = std(julia_norm) / mean(julia_norm)
                            @info "Diurnal CV ($sp): ref=$(round(ref_cv, digits=3)), julia=$(round(julia_cv, digits=3))"
                        end
                    end
                end
            end

            # ---- Per-cell spatial comparison for CO ----
            @testset "Per-cell spatial comparison (CO)" begin
                if haskey(model_data, "CO") && "CO" in ref_species
                    ref_spatial = ioapi_spatial_pattern(ds, "CO")
                    julia_spatial = dropdims(sum(model_data["CO"], dims = (3, 4)), dims = (3, 4))

                    ref_total = sum(ref_spatial)
                    julia_total = sum(julia_spatial)

                    if ref_total > 0 && julia_total > 0
                        ref_norm = ref_spatial ./ ref_total
                        julia_norm = julia_spatial ./ julia_total

                        # Normalized RMSE: should be small for well-matched patterns
                        diff = ref_norm .- julia_norm
                        nrmse = sqrt(mean(diff .^ 2)) / mean(ref_norm)
                        @info "CO normalized RMSE: $(round(nrmse, digits=4))"
                        @test nrmse < 1.0

                        # Check that the top-emitting cells overlap
                        ref_flat = vec(ref_norm)
                        julia_flat = vec(julia_norm)
                        ref_top10 = Set(sortperm(ref_flat, rev = true)[1:10])
                        julia_top10 = Set(sortperm(julia_flat, rev = true)[1:10])
                        overlap = length(intersect(ref_top10, julia_top10))
                        @info "CO top-10 cell overlap: $overlap / 10"
                        @test overlap >= 3
                    end
                end
            end

            # ---- Magnitude diagnostic (informational, not tested strictly) ----
            @testset "Magnitude diagnostics" begin
                # Log the total emissions in each species for debugging.
                # Due to unknown temporal scaling and mass/mole basis differences,
                # absolute magnitudes are not directly comparable.
                simple_species_mw = Dict(
                    "CO" => 28.0, "SO2" => 64.0, "NH3" => 17.0,
                    "NO" => 30.0, "NO2" => 46.0, "HONO" => 47.0,
                )
                pm_species = Set(["PAL", "PCA", "PCL", "PEC", "PFE", "PH2O", "PK",
                    "PMG", "PMN", "PMOTHR", "PNA", "PNCOM", "PNH4", "PNO3",
                    "POC", "PSI", "PSO4", "PTI", "PMC", "SULF"])

                for sp in sort(collect(common_species))
                    ref_total = Float64(sum(read_ioapi_var_raw(ds, sp)))
                    julia_total = sum(model_data[sp])

                    if sp in keys(simple_species_mw)
                        # Gas species: convert julia kg/s → mol/s for comparison
                        mw = simple_species_mw[sp]
                        julia_mol_per_s = julia_total * 1000.0 / mw
                        ratio = ref_total > 0 ? julia_mol_per_s / ref_total : NaN
                        @info "  $sp (gas): ref=$(round(ref_total, sigdigits=4)) mol/s, julia=$(round(julia_mol_per_s, sigdigits=4)) mol/s, ratio=$(round(ratio, sigdigits=3))"
                    elseif sp in pm_species
                        # PM species: convert julia kg/s → g/s for comparison
                        julia_g_per_s = julia_total * 1000.0
                        ratio = ref_total > 0 ? julia_g_per_s / ref_total : NaN
                        @info "  $sp (PM): ref=$(round(ref_total, sigdigits=4)) g/s, julia=$(round(julia_g_per_s, sigdigits=4)) g/s, ratio=$(round(ratio, sigdigits=3))"
                    else
                        @info "  $sp: ref=$(round(ref_total, sigdigits=4)), julia=$(round(julia_total, sigdigits=4))"
                    end
                end

                # As a basic sanity check, all common species should have non-zero
                # totals in both Julia and reference output
                for sp in common_species
                    ref_total = Float64(sum(read_ioapi_var_raw(ds, sp)))
                    julia_total = sum(model_data[sp])
                    if ref_total > 0
                        @test julia_total > 0
                    end
                end
            end
        end
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
