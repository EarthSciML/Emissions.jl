export NewPolygon, NewGridIrregular, setupSpatialProcessor, findCountyPolygon,
    GetIndex, recordToGrid, GridFactors, uniqueCoordinates, uniqueLoc,
    cell_bounds, cell_polygon, cell_area,
    build_regridder, grid_polygons

"""
    cell_bounds(grid::GridDef, idx::Int) -> (xmin, xmax, ymin, ymax)

Return the bounding box of grid cell `idx` as `(xmin, xmax, ymin, ymax)`.
"""
function cell_bounds(grid::GridDef, idx::Int)
    ext = grid.Extent[idx]
    xmin, ymin = ext[1]
    xmax, ymax = ext[2]
    return xmin, xmax, ymin, ymax
end

"""
    cell_polygon(grid::GridDef, idx::Int)

Construct a GeoInterface polygon for grid cell `idx` from its bounding box.
"""
function cell_polygon(grid::GridDef, idx::Int)
    xmin, xmax, ymin, ymax = cell_bounds(grid, idx)
    return GI.Polygon([[(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]])
end

"""
    cell_area(grid::GridDef, idx::Int) -> Float64

Return the area of grid cell `idx` computed from its bounding box.
"""
function cell_area(grid::GridDef, idx::Int)
    xmin, xmax, ymin, ymax = cell_bounds(grid, idx)
    return (xmax - xmin) * (ymax - ymin)
end

"""
    NewPolygon(coords::Vector{Tuple{Float64,Float64}})

Create a GeoInterface polygon from a vector of coordinate tuples.
The polygon is automatically closed if needed.
"""
function NewPolygon(coords::Vector{Tuple{Float64, Float64}})
    if coords[1] != coords[end]
        push!(coords, coords[1])
    end
    return GI.Polygon([coords])
end

"""
    NewGridIrregular(name, nx, ny, sr, dx, dy, x0, y0)

Create a `GridDef` with a regular (but potentially projected) grid.
"""
function NewGridIrregular(
        name::String, nx::Int, ny::Int, sr::String,
        dx::Float64, dy::Float64, x0::Float64, y0::Float64
    )
    extents = Vector{Tuple{Float64, Float64}}[]

    for j in 1:ny
        for i in 1:nx
            xmin = x0 + (i - 1) * dx
            xmax = x0 + i * dx
            ymin = y0 + (j - 1) * dy
            ymax = y0 + j * dy
            push!(extents, [(xmin, ymin), (xmax, ymax)])
        end
    end

    return GridDef(name, nx, ny, sr, extents)
end

"""
    NewGridIrregular(name, filepath, inputSR, outputSR)

Create an irregular grid by reading polygon coordinates from a file.
The file should contain polygon coordinate strings in the format used by the NEI processing notebook.
"""
function NewGridIrregular(name::String, filepath::String, inputSR::String, outputSR::String)
    extents = Vector{Tuple{Float64, Float64}}[]
    regex = r"X:(-?\d+(\.\d+)?(e[+\-]?\d+)?) Y:(-?\d+(\.\d+)?(e[+\-]?\d+)?)"
    trans = Proj.Transformation(inputSR, outputSR)

    open(filepath, "r") do file
        content = read(file, String)
        polygon_strs = split(replace(content, '\n' => ""), "]] [[")

        for polygon_str in polygon_strs
            cleaned_str = replace(polygon_str, r"[\[\]]" => "")
            matches = collect(eachmatch(regex, cleaned_str))

            if !isempty(matches)
                extent_x = Float64[]
                extent_y = Float64[]

                for match in matches
                    x1 = parse(Float64, match.captures[1])
                    y1 = parse(Float64, match.captures[4])
                    x, y = trans([x1, y1])
                    push!(extent_x, x)
                    push!(extent_y, y)
                end

                try
                    push!(
                        extents, [
                            (minimum(extent_x), minimum(extent_y)),
                            (maximum(extent_x), maximum(extent_y)),
                        ]
                    )
                catch e
                    @warn "Failed to create extent from coordinates: $e"
                    continue
                end
            end
        end
    end

    # Infer grid dimensions from number of cells
    ncells = length(extents)
    nx = isqrt(ncells)  # Approximate square grid
    ny = ncells ÷ nx + (ncells % nx > 0 ? 1 : 0)

    return GridDef(name, nx, ny, outputSR, extents)
end

"""
    setupSpatialProcessor(config::Config)

Set up a `SpatialProcessor` from a configuration object. Reads grid reference files,
surrogate specifications, and grid definitions.
"""
function setupSpatialProcessor(config::Config)
    gridRef = vcat([CSV.read(f, DataFrame; comment = "#") for f in config.f_gridRef]...)

    srgSpecs = readSrgSpecSMOKE(config.SrgSpec, config.SrgShapefileDirectory)

    # Create grid from file if it exists, otherwise create a default grid
    if isfile(config.GridFile)
        try
            grid = NewGridIrregular(config.GridName, config.GridFile, config.InputSR, config.OutputSR)
        catch e
            @warn "Failed to read grid file $(config.GridFile): $e. Using default 1x1 grid."
            grid = NewGridIrregular(config.GridName, 1, 1, config.OutputSR, 1.0, 1.0, 0.0, 0.0)
        end
    else
        @warn "Grid file $(config.GridFile) not found. Using default 1x1 grid."
        grid = NewGridIrregular(config.GridName, 1, 1, config.OutputSR, 1.0, 1.0, 0.0, 0.0)
    end

    return NewSpatialProcessor(srgSpecs, grid, gridRef, config.InputSR, false)
end

"""
    findCountyPolygon(fips::AbstractString, countyShapefile::AbstractString)

Find and return the polygon for a given FIPS code from a county shapefile.
Uses GeoDataFrames to read the shapefile. Returns a GeoInterface-compatible
geometry or `nothing` if not found.
"""
function findCountyPolygon(fips::AbstractString, countyShapefile::AbstractString)
    df = GeoDataFrames.read(countyShapefile)
    for row in eachrow(df)
        geofips = lpad(string(row.STATEFP) * string(row.COUNTYFP), 5, '0')
        if geofips == fips
            return row.geometry
        end
    end
    return nothing
end

"""
    GetIndex(lon, lat, grid::GridDef)

Find which grid cell contains the point (lon, lat) and return an `IndexInfo`.
Uses bounding box containment for efficient point-in-cell lookup.
"""
function GetIndex(lon::Float64, lat::Float64, grid::GridDef)
    rows = Int[]
    cols = Int[]
    fracs = Float64[]

    for j in 1:grid.Ny
        for i in 1:grid.Nx
            idx = (j - 1) * grid.Nx + i
            xmin, xmax, ymin, ymax = cell_bounds(grid, idx)
            if xmin <= lon <= xmax && ymin <= lat <= ymax
                push!(rows, j)
                push!(cols, i)
                push!(fracs, 1.0)
                return IndexInfo(rows, cols, fracs, true, true)
            end
        end
    end

    return IndexInfo(rows, cols, fracs, false, false)
end

"""
    GetIndex(geom, grid::GridDef)

Find which grid cells a geometry intersects, and return an `IndexInfo`
with the row/column indices and fractional coverage.
Uses `GeometryOps.intersection` to compute intersection polygons and areas
for each grid cell.
"""
function GetIndex(geom, grid::GridDef)
    rows = Int[]
    cols = Int[]
    fracs = Float64[]
    totalArea = GO.area(geom)

    if totalArea == 0.0
        return IndexInfo(rows, cols, fracs, false, false)
    end

    inGrid = false
    coveredArea = 0.0

    for j in 1:grid.Ny
        for i in 1:grid.Nx
            cell = cell_polygon(grid, (j - 1) * grid.Nx + i)
            inter = GO.intersection(geom, cell; target = GI.PolygonTrait())
            if !isempty(inter)
                area = sum(GO.area, inter)
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
            areas[j, i] = cell_area(grid, idx)
        end
    end
    return areas
end

"""
    uniqueCoordinates(lons::Vector{Float64}, lats::Vector{Float64})

Return the indices of unique (lon, lat) coordinate pairs.
"""
function uniqueCoordinates(lons::Vector{Float64}, lats::Vector{Float64})
    seen = Dict{Tuple{Float64, Float64}, Int}()
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
    seen = Dict{Tuple{Float64, Float64}, Int}()
    for i in eachindex(lons)
        key = (lons[i], lats[i])
        if !haskey(seen, key)
            seen[key] = i
        end
    end
    return seen
end

"""
    grid_polygons(grid::GridDef) -> Vector

Return a vector of GeoInterface polygons for all grid cells, ordered by
cell index `(j-1)*Nx + i`.
"""
function grid_polygons(grid::GridDef)
    polys = []
    for j in 1:grid.Ny
        for i in 1:grid.Nx
            idx = (j - 1) * grid.Nx + i
            push!(polys, cell_polygon(grid, idx))
        end
    end
    return polys
end

"""
    _to_vertex_ring(geom) -> Vector{Tuple{Float64,Float64}}

Extract the exterior ring vertex coordinates from a GeoInterface geometry
as a vector of (x,y) tuples suitable for ConservativeRegridding.
"""
function _to_vertex_ring(geom)
    ring = GI.getexterior(geom)
    coords = Tuple{Float64, Float64}[]
    for pt in GI.getpoint(ring)
        push!(coords, (Float64(GI.x(pt)), Float64(GI.y(pt))))
    end
    return coords
end

"""
    build_regridder(source_geometries::Vector, grid::GridDef) -> ConservativeRegridding.Regridder

Build a ConservativeRegridding.jl `Regridder` from source geometries to grid cells.

The regridder computes intersection areas between all source geometries and
grid cells using an efficient spatial index (STRtree). The resulting object
can be reused for multiple regridding operations on the same grid.

# Arguments
- `source_geometries`: Vector of GeoInterface-compatible polygon geometries.
- `grid::GridDef`: Target grid definition.

# Returns
A `ConservativeRegridding.Regridder` that can be used with `ConservativeRegridding.regrid`.
"""
function build_regridder(source_geometries::Vector, grid::GridDef)
    dst_verts = [_to_vertex_ring(p) for p in grid_polygons(grid)]
    src_verts = [_to_vertex_ring(g) for g in source_geometries]
    return ConservativeRegridding.Regridder(dst_verts, src_verts; normalize = false)
end
