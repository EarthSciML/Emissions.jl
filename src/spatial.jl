export NewPolygon, NewGridIrregular, setupSpatialProcessor, findCountyPolygon,
    GetIndex, recordToGrid, GridFactors, uniqueCoordinates, uniqueLoc

"""
    NewPolygon(coords::Vector{Tuple{Float64,Float64}})

Create a LibGEOS polygon from a vector of coordinate tuples.
The polygon is automatically closed if needed.
"""
function NewPolygon(coords::Vector{Tuple{Float64,Float64}})
    if coords[1] != coords[end]
        push!(coords, coords[1])
    end
    wkt = "POLYGON((" * join(["$(c[1]) $(c[2])" for c in coords], ",") * "))"
    return LibGEOS.readgeom(wkt)
end

"""
    NewGridIrregular(name, nx, ny, sr, dx, dy, x0, y0)

Create a `GridDef` with a regular (but potentially projected) grid.
"""
function NewGridIrregular(name::String, nx::Int, ny::Int, sr::String,
    dx::Float64, dy::Float64, x0::Float64, y0::Float64)
    cells = LibGEOS.Polygon[]
    extents = Vector{Tuple{Float64,Float64}}[]

    for j in 1:ny
        for i in 1:nx
            xmin = x0 + (i - 1) * dx
            xmax = x0 + i * dx
            ymin = y0 + (j - 1) * dy
            ymax = y0 + j * dy
            coords = [
                (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)
            ]
            wkt = "POLYGON((" * join(["$(c[1]) $(c[2])" for c in coords], ",") * "))"
            push!(cells, LibGEOS.readgeom(wkt))
            push!(extents, [(xmin, ymin), (xmax, ymax)])
        end
    end

    return GridDef(name, nx, ny, sr, cells, extents)
end

"""
    setupSpatialProcessor(config::Config)

Set up a `SpatialProcessor` from a configuration object. Reads grid reference files,
surrogate specifications, and grid definitions.
"""
function setupSpatialProcessor(config::Config)
    gridRef = vcat([CSV.read(f, DataFrame; comment="#") for f in config.f_gridRef]...)

    srgSpecs = readSrgSpecSMOKE(config.SrgSpec, config.SrgShapefileDirectory)

    # Read grid from file
    gridData = read_grid(config.GridFile)

    # Create a simple grid (placeholder - real implementation would parse grid file)
    grid = NewGridIrregular(config.GridName, 1, 1, config.OutputSR, 1.0, 1.0, 0.0, 0.0)

    return NewSpatialProcessor(srgSpecs, grid, gridRef, config.InputSR, false)
end

"""
    findCountyPolygon(fips::AbstractString, countyShapefile::AbstractString)

Find and return the polygon for a given FIPS code from a county shapefile.
"""
function findCountyPolygon(fips::AbstractString, countyShapefile::AbstractString)
    table = Shapefile.Table(countyShapefile)
    for row in table
        geofips = lpad(string(row.STATEFP) * string(row.COUNTYFP), 5, '0')
        if geofips == fips
            geom = Shapefile.shape(row)
            return geom
        end
    end
    return nothing
end

"""
    GetIndex(point_or_poly, grid::GridDef)

Find which grid cells a point or polygon intersects, and return an `IndexInfo`
with the row/column indices and fractional coverage.
"""
function GetIndex(lon::Float64, lat::Float64, grid::GridDef)
    point_wkt = "POINT($lon $lat)"
    point = LibGEOS.readgeom(point_wkt)

    rows = Int[]
    cols = Int[]
    fracs = Float64[]

    for j in 1:grid.Ny
        for i in 1:grid.Nx
            idx = (j - 1) * grid.Nx + i
            if LibGEOS.contains(grid.Cells[idx], point)
                push!(rows, j)
                push!(cols, i)
                push!(fracs, 1.0)
                return IndexInfo(rows, cols, fracs, true, true)
            end
        end
    end

    return IndexInfo(rows, cols, fracs, false, false)
end

function GetIndex(poly::LibGEOS.Polygon, grid::GridDef)
    rows = Int[]
    cols = Int[]
    fracs = Float64[]
    totalArea = LibGEOS.area(poly)

    if totalArea == 0.0
        return IndexInfo(rows, cols, fracs, false, false)
    end

    inGrid = false
    coveredArea = 0.0

    for j in 1:grid.Ny
        for i in 1:grid.Nx
            idx = (j - 1) * grid.Nx + i
            if LibGEOS.intersects(grid.Cells[idx], poly)
                intersection = LibGEOS.intersection(grid.Cells[idx], poly)
                area = LibGEOS.area(intersection)
                if area > 0.0
                    push!(rows, j)
                    push!(cols, i)
                    push!(fracs, area / totalArea)
                    coveredArea += area
                    inGrid = true
                end
            end
        end
    end

    coveredByGrid = coveredArea / totalArea > 0.99

    return IndexInfo(rows, cols, fracs, inGrid, coveredByGrid)
end

"""
    recordToGrid(emissions::Float64, index::IndexInfo, grid_nx::Int, grid_ny::Int)

Distribute emissions to grid cells according to the fractional coverage in `IndexInfo`.
Returns a sparse matrix of emissions.
"""
function recordToGrid(emissions::Float64, index::IndexInfo, grid_nx::Int, grid_ny::Int)
    result = spzeros(grid_ny, grid_nx)
    for k in eachindex(index.rows)
        result[index.rows[k], index.cols[k]] += emissions * index.fracs[k]
    end
    return result
end

"""
    GridFactors(grid::GridDef)

Return a matrix of grid cell areas.
"""
function GridFactors(grid::GridDef)
    areas = zeros(grid.Ny, grid.Nx)
    for j in 1:grid.Ny
        for i in 1:grid.Nx
            idx = (j - 1) * grid.Nx + i
            areas[j, i] = LibGEOS.area(grid.Cells[idx])
        end
    end
    return areas
end

"""
    uniqueCoordinates(lons::Vector{Float64}, lats::Vector{Float64})

Return the indices of unique (lon, lat) coordinate pairs.
"""
function uniqueCoordinates(lons::Vector{Float64}, lats::Vector{Float64})
    seen = Dict{Tuple{Float64,Float64},Int}()
    indices = Int[]
    for i in eachindex(lons)
        key = (lons[i], lats[i])
        if !haskey(seen, key)
            seen[key] = i
            push!(indices, i)
        end
    end
    return indices
end

"""
    uniqueLoc(lons::Vector{Float64}, lats::Vector{Float64})

Return a dictionary mapping unique (lon, lat) pairs to their first occurrence index.
"""
function uniqueLoc(lons::Vector{Float64}, lats::Vector{Float64})
    seen = Dict{Tuple{Float64,Float64},Int}()
    for i in eachindex(lons)
        key = (lons[i], lats[i])
        if !haskey(seen, key)
            seen[key] = i
        end
    end
    return seen
end
