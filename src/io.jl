export strip_missing, getCountry, read_grid, getShapefilePath,
    validateShapefile, readSrgSpecSMOKE, NewSpatialProcessor

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
    read_grid(file::AbstractString)

Read a grid definition file (IOAPI-style GRIDDESC format) and return a DataFrame
with columns FIPS, Longitude, and Latitude.
"""
function read_grid(file::AbstractString)
    lines = readlines(file)
    records = DataFrame(FIPS=String[], Longitude=Float64[], Latitude=Float64[])
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
            push!(records, (FIPS=fips, Longitude=lon, Latitude=lat))
        end
    end
    return records
end

"""
    getShapefilePath(dir::AbstractString, name::AbstractString, check::Bool=true)

Find the full path to a shapefile within the given directory tree.
Searches recursively for files with `.shp` extension matching the given name.
"""
function getShapefilePath(dir::AbstractString, name::AbstractString, check::Bool=true)
    if endswith(name, ".shp")
        base = name[1:end-4]
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
function readSrgSpecSMOKE(file::AbstractString, SrgShapefileDirectory::AbstractString, CheckShapefiles::Bool=false)
    df = CSV.read(file, DataFrame; header=false, comment="#", silencewarnings=true)
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

        push!(srgs, SurrogateSpec(Region, Name, Code, DataShapefile, DataAttribute,
            WeightShapefile, Details, BackupSurrogateNames,
            WeightColumns, WeightFactors, FilterFunction,
            MergeNames, MergeMultipliers))
    end

    return srgs
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
        100,
        10,
    )
    return sp
end
