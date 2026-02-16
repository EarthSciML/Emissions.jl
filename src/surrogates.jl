export generate_data_sparse_matrices, generate_weight_sparse_matrices,
    generate_grid_sparse_matrices, generate_countySurrogate, update_locIndex,
    generate_data_sparse_matrices_cr, generate_weight_sparse_matrices_cr

"""
    find_column_name(df::DataFrame, name::AbstractString)

Find a column name in a DataFrame, case-insensitive.
Returns the actual column name or throws an error.
"""
function find_column_name(df::DataFrame, name::AbstractString)
    for col in names(df)
        if lowercase(col) == lowercase(name)
            return col
        end
    end
    error("Column '$name' not found in DataFrame")
end

"""
    read_crs_epsg(shapefile_path::AbstractString)

Read the CRS from a shapefile's .prj file and return the EPSG code as an integer.
Uses Proj.jl for CRS identification.
"""
function read_crs_epsg(shapefile_path::AbstractString)
    prj_path = replace(shapefile_path, r"\.shp$" => ".prj")
    if !isfile(prj_path)
        error("No .prj file found for shapefile: $shapefile_path")
    end
    wkt = read(prj_path, String)
    crs = Proj.CRS(wkt)
    identified = Proj.identify(crs; auth_name = "EPSG")
    if isempty(identified)
        error("Could not identify EPSG code for CRS in $prj_path")
    end
    # Use the crs field from the NamedTuple returned by identify
    crs_obj = identified[1].crs
    code_str = Proj.proj_get_id_code(crs_obj.pj, 0)
    return parse(Int, code_str)
end

"""
    generate_data_sparse_matrices(shapefile_path, attribute, grid, inputSR)

Read a data shapefile using GeoDataFrames and generate sparse matrices representing the
spatial allocation of each region (identified by `attribute`) to grid cells.
Returns a dictionary mapping attribute values to sparse matrices.
"""
function generate_data_sparse_matrices(
        shapefile_path::AbstractString,
        attribute::AbstractString, grid::GridDef, inputSR::String
    )

    df = GeoDataFrames.read(shapefile_path)
    col = find_column_name(df, attribute)

    result = Dict{String, SparseMatrixCSC{Float64, Int}}()

    for row in eachrow(df)
        key = string(row[Symbol(col)])
        geom = row.geometry

        for j in 1:grid.Ny
            for i in 1:grid.Nx
                cell = cell_polygon(grid, (j - 1) * grid.Nx + i)
                inter = GO.intersection(geom, cell; target = GI.PolygonTrait())
                if !isempty(inter)
                    area = sum(GO.area, inter)
                    if area > 0.0
                        if !haskey(result, key)
                            result[key] = spzeros(grid.Ny, grid.Nx)
                        end
                        result[key][j, i] += area
                    end
                end
            end
        end
    end

    return result
end

"""
    generate_weight_sparse_matrices(shapefile_path, weight_columns, weight_factors, grid, inputSR)

Read a weight shapefile using GeoDataFrames and generate sparse matrices
of weighted values on the grid.
"""
function generate_weight_sparse_matrices(
        shapefile_path::AbstractString,
        weight_columns::Vector{String}, weight_factors::Vector{Float64},
        grid::GridDef, inputSR::String
    )

    df = GeoDataFrames.read(shapefile_path)

    result = spzeros(grid.Ny, grid.Nx)

    for row in eachrow(df)
        weight = 0.0
        for (col, factor) in zip(weight_columns, weight_factors)
            actual_col = find_column_name(df, col)
            val = row[Symbol(actual_col)]
            if !ismissing(val) && val isa Number
                weight += val * factor
            end
        end

        geom = row.geometry
        for j in 1:grid.Ny
            for i in 1:grid.Nx
                cell = cell_polygon(grid, (j - 1) * grid.Nx + i)
                inter = GO.intersection(geom, cell; target = GI.PolygonTrait())
                if !isempty(inter)
                    area = sum(GO.area, inter)
                    if area > 0.0
                        result[j, i] += weight * area
                    end
                end
            end
        end
    end

    return result
end

"""
    generate_grid_sparse_matrices(grid::GridDef)

Generate a sparse matrix of grid cell areas.
"""
function generate_grid_sparse_matrices(grid::GridDef)
    result = spzeros(grid.Ny, grid.Nx)
    for j in 1:grid.Ny
        for i in 1:grid.Nx
            idx = (j - 1) * grid.Nx + i
            result[j, i] = cell_area(grid, idx)
        end
    end
    return result
end

"""
    generate_countySurrogate(data_matrices, weight_matrix, grid_matrix)

Generate a surrogate matrix by combining data and weight matrices.
For each region in data_matrices, the surrogate is computed as:
  surrogate = (data .* weight) / sum(data .* weight)
Falls back to area-based allocation if weight sum is zero.
"""
function generate_countySurrogate(
        data_matrices::Dict{String, SparseMatrixCSC{Float64, Int}},
        weight_matrix::SparseMatrixCSC{Float64, Int},
        grid_matrix::SparseMatrixCSC{Float64, Int}
    )

    result = Dict{String, SparseMatrixCSC{Float64, Int}}()

    for (key, data) in data_matrices
        weighted = data .* weight_matrix
        total = sum(weighted)

        if total > 0.0
            result[key] = weighted ./ total
        else
            # Fall back to area-based allocation
            area_total = sum(data)
            if area_total > 0.0
                result[key] = data ./ area_total
            else
                result[key] = spzeros(size(data)...)
            end
        end
    end

    return result
end

"""
    update_locIndex(locIndex::Dict, fips::AbstractString, surrogates::Dict)

Update a location index dictionary with surrogate allocation data for a given FIPS code.
Returns the IndexInfo for the FIPS code from the surrogates, or a default empty IndexInfo.
"""
function update_locIndex(
        locIndex::Dict{String, IndexInfo}, fips::AbstractString,
        surrogates::Dict{String, SparseMatrixCSC{Float64, Int}}
    )

    if haskey(surrogates, fips)
        srg = surrogates[fips]
        rows_idx, cols_idx, vals = findnz(srg)
        idx = IndexInfo(rows_idx, cols_idx, vals, !isempty(rows_idx), !isempty(rows_idx))
        locIndex[fips] = idx
        return idx
    end

    # Return empty index
    idx = IndexInfo(Int[], Int[], Float64[], false, false)
    locIndex[fips] = idx
    return idx
end

"""
    generate_data_sparse_matrices_cr(shapefile_path, attribute, grid, inputSR)

Generate data sparse matrices using ConservativeRegridding.jl for efficient
area-overlap computation via spatial indexing.

Same interface as [`generate_data_sparse_matrices`](@ref), but uses
precomputed intersection area matrices for better performance with large grids.
"""
function generate_data_sparse_matrices_cr(
        shapefile_path::AbstractString,
        attribute::AbstractString, grid::GridDef, inputSR::String
    )
    df = GeoDataFrames.read(shapefile_path)
    col = find_column_name(df, attribute)

    # Group geometries by attribute value
    groups = Dict{String, Vector{Int}}()
    geoms = []
    for (i, row) in enumerate(eachrow(df))
        key = string(row[Symbol(col)])
        push!(geoms, row.geometry)
        if !haskey(groups, key)
            groups[key] = Int[]
        end
        push!(groups[key], i)
    end

    if isempty(geoms)
        return Dict{String, SparseMatrixCSC{Float64, Int}}()
    end

    # Build regridder: grid cells (dst) × source polygons (src)
    dst_verts = [_to_vertex_ring(p) for p in grid_polygons(grid)]
    src_verts = [_to_vertex_ring(g) for g in geoms]
    regridder = ConservativeRegridding.Regridder(dst_verts, src_verts; normalize = false)
    areas = regridder.intersections  # dst_cells × src_polygons

    result = Dict{String, SparseMatrixCSC{Float64, Int}}()
    ncells = grid.Ny * grid.Nx

    for (key, src_indices) in groups
        mat = spzeros(grid.Ny, grid.Nx)
        for src_idx in src_indices
            for cell_idx in 1:ncells
                a = areas[cell_idx, src_idx]
                if a > 0.0
                    j = (cell_idx - 1) ÷ grid.Nx + 1
                    i = (cell_idx - 1) % grid.Nx + 1
                    mat[j, i] += a
                end
            end
        end
        result[key] = mat
    end

    return result
end

"""
    generate_weight_sparse_matrices_cr(shapefile_path, weight_columns, weight_factors, grid, inputSR)

Generate weight sparse matrices using ConservativeRegridding.jl for efficient
area-overlap computation.

Same interface as [`generate_weight_sparse_matrices`](@ref), but uses
precomputed intersection area matrices for better performance.
"""
function generate_weight_sparse_matrices_cr(
        shapefile_path::AbstractString,
        weight_columns::Vector{String}, weight_factors::Vector{Float64},
        grid::GridDef, inputSR::String
    )
    df = GeoDataFrames.read(shapefile_path)

    # Compute weights per source polygon
    n_src = nrow(df)
    weights = zeros(Float64, n_src)
    geoms = []

    for (i, row) in enumerate(eachrow(df))
        w = 0.0
        for (col, factor) in zip(weight_columns, weight_factors)
            actual_col = find_column_name(df, col)
            val = row[Symbol(actual_col)]
            if !ismissing(val) && val isa Number
                w += val * factor
            end
        end
        weights[i] = w
        push!(geoms, row.geometry)
    end

    if isempty(geoms)
        return spzeros(grid.Ny, grid.Nx)
    end

    # Build regridder
    dst_verts = [_to_vertex_ring(p) for p in grid_polygons(grid)]
    src_verts = [_to_vertex_ring(g) for g in geoms]
    regridder = ConservativeRegridding.Regridder(dst_verts, src_verts; normalize = false)
    areas = regridder.intersections  # dst_cells × src_polygons

    result = spzeros(grid.Ny, grid.Nx)
    ncells = grid.Ny * grid.Nx

    for src_idx in 1:n_src
        w = weights[src_idx]
        w == 0.0 && continue
        for cell_idx in 1:ncells
            a = areas[cell_idx, src_idx]
            if a > 0.0
                j = (cell_idx - 1) ÷ grid.Nx + 1
                i = (cell_idx - 1) % grid.Nx + 1
                result[j, i] += w * a
            end
        end
    end

    return result
end
