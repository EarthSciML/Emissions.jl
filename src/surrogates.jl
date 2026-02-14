export generate_data_sparse_matrices, generate_weight_sparse_matrices,
    generate_grid_sparse_matrices, generate_countySurrogate, update_locIndex

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
Uses Proj.jl for CRS identification instead of Python pyproj.
"""
function read_crs_epsg(shapefile_path::AbstractString)
    prj_path = replace(shapefile_path, r"\.shp$" => ".prj")
    if !isfile(prj_path)
        error("No .prj file found for shapefile: $shapefile_path")
    end
    wkt = read(prj_path, String)
    crs = Proj.CRS(wkt)
    identified = Proj.identify(crs; auth_name="EPSG")
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

Read a data shapefile and generate sparse matrices representing the
spatial allocation of each region (identified by `attribute`) to grid cells.
Returns a dictionary mapping attribute values to sparse matrices.
"""
function generate_data_sparse_matrices(shapefile_path::AbstractString,
    attribute::AbstractString, grid::GridDef, inputSR::String)

    table = Shapefile.Table(shapefile_path)
    col = find_column_name(DataFrame(table), attribute)

    result = Dict{String,SparseMatrixCSC{Float64,Int}}()

    for row in table
        key = string(getproperty(row, Symbol(col)))
        geom = Shapefile.shape(row)

        # Convert geometry to LibGEOS
        wkt_str = GeoInterface.astext(geom)
        lgeos_geom = LibGEOS.readgeom(wkt_str)

        for j in 1:grid.Ny
            for i in 1:grid.Nx
                idx = (j - 1) * grid.Nx + i
                if LibGEOS.intersects(grid.Cells[idx], lgeos_geom)
                    intersection = LibGEOS.intersection(grid.Cells[idx], lgeos_geom)
                    area = LibGEOS.area(intersection)
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

Read a weight shapefile and generate sparse matrices of weighted values on the grid.
"""
function generate_weight_sparse_matrices(shapefile_path::AbstractString,
    weight_columns::Vector{String}, weight_factors::Vector{Float64},
    grid::GridDef, inputSR::String)

    table = Shapefile.Table(shapefile_path)
    df = DataFrame(table)

    result = spzeros(grid.Ny, grid.Nx)

    for row in table
        geom = Shapefile.shape(row)
        wkt_str = GeoInterface.astext(geom)
        lgeos_geom = LibGEOS.readgeom(wkt_str)

        weight = 0.0
        for (col, factor) in zip(weight_columns, weight_factors)
            actual_col = find_column_name(df, col)
            val = getproperty(row, Symbol(actual_col))
            if !ismissing(val) && val isa Number
                weight += val * factor
            end
        end

        for j in 1:grid.Ny
            for i in 1:grid.Nx
                idx = (j - 1) * grid.Nx + i
                if LibGEOS.intersects(grid.Cells[idx], lgeos_geom)
                    intersection = LibGEOS.intersection(grid.Cells[idx], lgeos_geom)
                    area = LibGEOS.area(intersection)
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
            result[j, i] = LibGEOS.area(grid.Cells[idx])
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
function generate_countySurrogate(data_matrices::Dict{String,SparseMatrixCSC{Float64,Int}},
    weight_matrix::SparseMatrixCSC{Float64,Int},
    grid_matrix::SparseMatrixCSC{Float64,Int})

    result = Dict{String,SparseMatrixCSC{Float64,Int}}()

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
function update_locIndex(locIndex::Dict{String,IndexInfo}, fips::AbstractString,
    surrogates::Dict{String,SparseMatrixCSC{Float64,Int}})

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
