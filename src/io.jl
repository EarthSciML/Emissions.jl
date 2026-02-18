export strip_missing, getCountry, normalize_country, read_grid, read_gridref,
    getShapefilePath, validateShapefile, readSrgSpecSMOKE, NewSpatialProcessor,
    read_griddesc

"""
    strip_missing(df::DataFrame)

Replace all `missing` values in a DataFrame with empty strings.
"""
function strip_missing(df::DataFrame)
    for col in names(df)
        df[!, col] = [ismissing(v) ? "" : v for v in df[!, col]]
    end
    return df
end

"""
    getCountry(fips::AbstractString)

Return the country name based on the first digit of the FIPS code.
- '0' or '9' => "US"
- '1' => "Mexico"
- '2' => "Canada"
"""
function getCountry(fips::AbstractString)
    if length(fips) == 0
        return "Unknown"
    end
    c = fips[1]
    if c == '0' || c == '9'
        return "US"
    elseif c == '1'
        return "Mexico"
    elseif c == '2'
        return "Canada"
    else
        return "Unknown"
    end
end

"""
    normalize_country(country::AbstractString) -> String

Normalize country codes to standard three-letter format.
Maps "US" -> "USA", "0" -> "USA", "1" -> "Canada", "2" -> "Mexico".
Already-standard codes are passed through unchanged.
"""
function normalize_country(country::AbstractString)
    c = strip(country)
    if c == "US" || c == "0"
        return "USA"
    elseif c == "1"
        return "Canada"
    elseif c == "2"
        return "Mexico"
    else
        return String(c)
    end
end

"""
    read_gridref(file_path::AbstractString) -> DataFrame

Read a SMOKE-format semicolon-delimited grid reference file.
Lines starting with `#` are skipped. Each data line has the format:
`FIPS;SCC;Surrogate!comment`

Extracts COUNTRY from 6-digit FIPS codes and normalizes via [`normalize_country`](@ref).
Returns a DataFrame with columns: `[:COUNTRY, :FIPS, :SCC, :Surrogate]`.
"""
function read_gridref(file_path::AbstractString)
    lines = readlines(file_path)
    records = DataFrame(COUNTRY = String[], FIPS = String[], SCC = String[], Surrogate = Int[])
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        # Remove inline comments after '!'
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ';')
        if length(parts) >= 3
            fips_raw = strip(parts[1])
            scc = strip(parts[2])
            surrogate = parse(Int, strip(parts[3]))

            # Extract country from 6-digit FIPS
            if length(fips_raw) == 6
                country_digit = string(fips_raw[1])
                fips = lpad(fips_raw[2:end], 5, '0')
            else
                country_digit = "0"
                fips = lpad(fips_raw, 5, '0')
            end
            country = normalize_country(country_digit)
            push!(records, (COUNTRY = country, FIPS = fips, SCC = scc, Surrogate = surrogate))
        end
    end
    return records
end

"""
    read_grid(file::AbstractString)

Read a CSV grid definition file with columns for FIPS, Longitude, and Latitude.
Returns a DataFrame with columns `FIPS`, `Longitude`, and `Latitude`.

For IOAPI-format GRIDDESC files, use [`read_griddesc`](@ref) instead.
"""
function read_grid(file::AbstractString)
    lines = readlines(file)
    records = DataFrame(FIPS = String[], Longitude = Float64[], Latitude = Float64[])
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        parts = split(line, ',')
        if length(parts) >= 3
            fips = strip(parts[1])
            # Strip leading country digit if 6-digit FIPS
            if length(fips) == 6
                fips = fips[2:end]
            end
            fips = lpad(fips, 5, '0')
            lon = parse(Float64, strip(parts[2]))
            lat = parse(Float64, strip(parts[3]))
            push!(records, (FIPS = fips, Longitude = lon, Latitude = lat))
        end
    end
    return records
end

"""
    getShapefilePath(dir::AbstractString, name::AbstractString, check::Bool=true)

Find the full path to a shapefile within the given directory tree.
Searches recursively for files with `.shp` extension matching the given name.
"""
function getShapefilePath(dir::AbstractString, name::AbstractString, check::Bool = true)
    if endswith(name, ".shp")
        base = name[1:(end - 4)]
    else
        base = name
    end

    for (root, dirs, files) in walkdir(dir)
        for f in files
            if endswith(f, ".shp") && startswith(f, base)
                return joinpath(root, f)
            end
        end
    end

    if check
        @warn "Shapefile not found: $name in $dir"
    end
    return ""
end

"""
    validateShapefile(path::AbstractString)

Check that a shapefile path is non-empty and the file exists.
"""
function validateShapefile(path::AbstractString)
    return !isempty(path) && isfile(path)
end

"""
    readSrgSpecSMOKE(file::AbstractString, SrgShapefileDirectory::AbstractString, CheckShapefiles::Bool=false)

Read a SMOKE-format spatial surrogate specification file and return
a vector of `SurrogateSpec` objects.
"""
function readSrgSpecSMOKE(file::AbstractString, SrgShapefileDirectory::AbstractString, CheckShapefiles::Bool = false)
    df = CSV.read(file, DataFrame; header = false, comment = "#", silencewarnings = true)
    strip_missing(df)

    srgs = SurrogateSpec[]

    for row in eachrow(df)
        vals = [string(row[c]) for c in names(df)]

        Region = strip(get(vals, 1, ""))
        Name = strip(get(vals, 2, ""))
        Code = parse(Int, strip(get(vals, 3, "0")))
        DataShapefile = strip(get(vals, 4, ""))
        DataAttribute = strip(get(vals, 5, ""))
        WeightShapefile = strip(get(vals, 6, ""))
        Details = strip(get(vals, 7, ""))

        BackupSurrogateNames = String[]
        backups = strip(get(vals, 8, ""))
        if !isempty(backups)
            BackupSurrogateNames = [strip(s) for s in split(backups, ";") if !isempty(strip(s))]
        end

        WeightColumns = String[]
        wcols = strip(get(vals, 9, ""))
        if !isempty(wcols)
            WeightColumns = [strip(s) for s in split(wcols, ";") if !isempty(strip(s))]
        end

        WeightFactors = Float64[]
        wfactors = strip(get(vals, 10, ""))
        if !isempty(wfactors)
            WeightFactors = [parse(Float64, strip(s)) for s in split(wfactors, ";") if !isempty(strip(s))]
        end

        FilterFunction = strip(get(vals, 11, ""))

        MergeNames = String[]
        mnames = strip(get(vals, 12, ""))
        if !isempty(mnames)
            MergeNames = [strip(s) for s in split(mnames, ";") if !isempty(strip(s))]
        end

        MergeMultipliers = Float64[]
        mmults = strip(get(vals, 13, ""))
        if !isempty(mmults)
            MergeMultipliers = [parse(Float64, strip(s)) for s in split(mmults, ";") if !isempty(strip(s))]
        end

        if CheckShapefiles && !isempty(DataShapefile)
            DataShapefile = getShapefilePath(SrgShapefileDirectory, String(DataShapefile), CheckShapefiles)
            validateShapefile(DataShapefile)
        end

        if !isempty(WeightShapefile)
            WeightShapefile = getShapefilePath(SrgShapefileDirectory, String(WeightShapefile), CheckShapefiles)
            isValidWeightShapefile = validateShapefile(WeightShapefile)
        end

        push!(
            srgs, SurrogateSpec(
                Region, Name, Code, DataShapefile, DataAttribute,
                WeightShapefile, Details, BackupSurrogateNames,
                WeightColumns, WeightFactors, FilterFunction,
                MergeNames, MergeMultipliers
            )
        )
    end

    return srgs
end

"""
    _coordtype_to_proj4(coordtype, P_ALP, P_BET, P_GAM, XCENT, YCENT) -> String

Convert an IOAPI coordinate type code and projection parameters to a Proj4 string.

Supported coordinate types:
- 1 (LATGRD3): Geographic latitude-longitude
- 2 (LAMGRD3): Lambert Conformal Conic
- 5 (UTMGRD3): Universal Transverse Mercator
- 6 (POLGRD3): Polar Stereographic
- 7 (EQMGRD3): Equatorial Mercator
- 8 (TRMGRD3): Transverse Mercator
"""
function _coordtype_to_proj4(
        coordtype::Int, P_ALP::Float64, P_BET::Float64,
        P_GAM::Float64, XCENT::Float64, YCENT::Float64
    )
    if coordtype == 1
        return "+proj=longlat +datum=WGS84"
    elseif coordtype == 2
        return "+proj=lcc +lat_1=$(P_ALP) +lat_2=$(P_BET) +lon_0=$(P_GAM) +lat_0=$(YCENT) +datum=WGS84"
    elseif coordtype == 5
        zone = round(Int, P_ALP)
        return "+proj=utm +zone=$(zone) +datum=WGS84"
    elseif coordtype == 6
        return "+proj=stere +lat_0=$(P_ALP) +lon_0=$(P_GAM) +datum=WGS84"
    elseif coordtype == 7
        return "+proj=merc +lon_0=$(P_GAM) +datum=WGS84"
    elseif coordtype == 8
        return "+proj=tmerc +lon_0=$(P_GAM) +lat_0=$(YCENT) +datum=WGS84"
    else
        error("Unsupported IOAPI coordinate type: $coordtype")
    end
end

"""
    read_griddesc(filepath::AbstractString, grid_name::AbstractString) -> GridDef

Read an IOAPI-format GRIDDESC file and return a `GridDef` for the named grid.

The GRIDDESC format has two segments separated by blank lines or a line containing
only whitespace:

**Segment 1** — Coordinate system definitions (between first and second `' '` delimiters):
Each entry is two lines:
1. Coordinate system name (quoted string)
2. Seven numeric fields: `coordtype P_ALP P_BET P_GAM XCENT YCENT`

**Segment 2** — Grid definitions (between second and third `' '` delimiters):
Each entry is two lines:
1. Grid name (quoted string)
2. Nine fields: `coord_name XORIG YORIG XCELL YCELL NCOLS NROWS NTHIK`

The function looks up the requested `grid_name`, resolves its coordinate system,
converts the IOAPI coordinate type to a Proj4 string, and constructs a `GridDef`
via [`NewGridIrregular`](@ref).

# Arguments
- `filepath`: Path to the GRIDDESC file.
- `grid_name`: Name of the grid to extract (must match a grid name in Segment 2).

# Returns
A `GridDef` object for the specified grid.

# Example
```julia
grid = read_griddesc("GRIDDESC", "US36_CRO")
```
"""
function read_griddesc(filepath::AbstractString, grid_name::AbstractString)
    lines = readlines(filepath)

    # Parse into non-empty, non-comment lines
    cleaned = String[]
    for line in lines
        s = strip(line)
        if !isempty(s) && !startswith(s, "#")
            push!(cleaned, s)
        end
    end

    # Split into segments by delimiter lines (lines that are just a quoted space or similar sentinel)
    # GRIDDESC uses lines containing only ' ' (with possible whitespace) as segment delimiters
    segments = Vector{Vector{String}}()
    current_segment = String[]
    for line in cleaned
        # A delimiter line is one that contains only quotes and whitespace, e.g., "' '" or "''"
        stripped = replace(line, r"['\"\s]" => "")
        if isempty(stripped)
            if !isempty(current_segment)
                push!(segments, current_segment)
                current_segment = String[]
            end
        else
            push!(current_segment, line)
        end
    end
    if !isempty(current_segment)
        push!(segments, current_segment)
    end

    if length(segments) < 2
        error("GRIDDESC file must contain at least two segments (coordinate systems and grids)")
    end

    # Parse Segment 1: Coordinate systems
    # Each entry is a pair of lines: name line, then parameters line
    coord_systems = Dict{
        String, NamedTuple{
            (:coordtype, :P_ALP, :P_BET, :P_GAM, :XCENT, :YCENT),
            Tuple{Int, Float64, Float64, Float64, Float64, Float64},
        },
    }()

    seg1 = segments[1]
    i = 1
    while i + 1 <= length(seg1)
        coord_name = replace(strip(seg1[i]), r"['\"]" => "")
        params = split(strip(seg1[i + 1]))
        if length(params) >= 6
            coordtype = parse(Int, params[1])
            P_ALP = parse(Float64, params[2])
            P_BET = parse(Float64, params[3])
            P_GAM = parse(Float64, params[4])
            XCENT = parse(Float64, params[5])
            YCENT = parse(Float64, params[6])
            coord_systems[coord_name] = (
                coordtype = coordtype, P_ALP = P_ALP, P_BET = P_BET,
                P_GAM = P_GAM, XCENT = XCENT, YCENT = YCENT,
            )
        end
        i += 2
    end

    # Parse Segment 2: Grid definitions
    seg2 = segments[2]
    i = 1
    while i + 1 <= length(seg2)
        gname = replace(strip(seg2[i]), r"['\"]" => "")
        params = split(strip(seg2[i + 1]))

        if gname == grid_name && length(params) >= 8
            coord_name = replace(strip(params[1]), r"['\"]" => "")
            XORIG = parse(Float64, params[2])
            YORIG = parse(Float64, params[3])
            XCELL = parse(Float64, params[4])
            YCELL = parse(Float64, params[5])
            NCOLS = parse(Int, params[6])
            NROWS = parse(Int, params[7])
            # NTHIK = parse(Int, params[8])  # boundary thickness, not needed for GridDef

            if !haskey(coord_systems, coord_name)
                error("Coordinate system '$coord_name' referenced by grid '$grid_name' not found in GRIDDESC file")
            end

            cs = coord_systems[coord_name]
            proj4 = _coordtype_to_proj4(
                cs.coordtype, cs.P_ALP, cs.P_BET,
                cs.P_GAM, cs.XCENT, cs.YCENT
            )

            return NewGridIrregular(
                grid_name, NCOLS, NROWS, proj4,
                XCELL, YCELL, XORIG, YORIG
            )
        end
        i += 2
    end

    error("Grid '$grid_name' not found in GRIDDESC file: $filepath")
end

"""
    NewSpatialProcessor(srgSpecs, grids, gridRef, inputSR, matchFullSCC)

Create a new `SpatialProcessor` with default cache and merge depth settings.
"""
function NewSpatialProcessor(srgSpecs::Vector{SurrogateSpec}, grids::GridDef, gridRef::DataFrame, inputSR::String, matchFullSCC::Bool)
    sp = SpatialProcessor(
        srgSpecs,
        grids,
        gridRef,
        inputSR,
        matchFullSCC,
        10,
    )
    return sp
end
