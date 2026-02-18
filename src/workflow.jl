using Printf: @sprintf

export location_key, compute_grid_indices, refine_indices_with_surrogates,
    allocate_emissions_to_grid, process_emissions_spatial

"""
    location_key(fips::AbstractString, lon, lat) -> String

Compute a unique location key for an emissions record.

Point sources with valid coordinates return `"FIPS_lon_lat"`.
Area sources (missing coordinates) return the FIPS code.

This key is used to look up pre-computed grid cell allocations in the
location index dictionary returned by [`compute_grid_indices`](@ref).

# Examples
```jldoctest
julia> location_key("36001", -73.5, 40.5)
"36001_-73.500000_40.500000"

julia> location_key("36001", missing, missing)
"36001"
```
"""
function location_key(fips::AbstractString, lon, lat)
    if !ismissing(lon) && !ismissing(lat)
        return @sprintf("%s_%.6f_%.6f", fips, Float64(lon), Float64(lat))
    end
    return String(fips)
end

"""
    compute_grid_indices(emissions::DataFrame, grid::GridDef;
                         counties_shapefile::String="") -> Dict{String, IndexInfo}

Compute grid cell indices for all unique emission locations.

For each unique location in the emissions data:
- **Point sources** (rows with non-missing `LONGITUDE` and `LATITUDE`): finds the
  grid cell containing the point using [`GetIndex`](@ref).
- **Area sources** (without coordinates): uses county polygon intersection with
  the grid if `counties_shapefile` is provided. The county polygon is looked up
  by FIPS code using [`findCountyPolygon`](@ref).

Returns a `Dict{String, IndexInfo}` mapping location keys (see [`location_key`](@ref))
to their grid cell allocations.

# Arguments
- `emissions::DataFrame`: Must have a `:FIPS` column. May have `:LONGITUDE` and `:LATITUDE`.
- `grid::GridDef`: Target grid definition.
- `counties_shapefile::String=""`: Path to counties shapefile for area source indexing.
"""
function compute_grid_indices(
        emissions::DataFrame, grid::GridDef;
        counties_shapefile::String = ""
    )
    locIndex = Dict{String, IndexInfo}()
    has_lon = hasproperty(emissions, :LONGITUDE)
    has_lat = hasproperty(emissions, :LATITUDE)

    for row in eachrow(emissions)
        fips = string(row.FIPS)
        lon = has_lon ? row.LONGITUDE : missing
        lat = has_lat ? row.LATITUDE : missing
        key = location_key(fips, lon, lat)

        haskey(locIndex, key) && continue

        is_point = !ismissing(lon) && !ismissing(lat)

        if is_point
            idx = GetIndex(Float64(lon), Float64(lat), grid)
            locIndex[key] = idx
        elseif !isempty(counties_shapefile) && isfile(counties_shapefile)
            county_geom = findCountyPolygon(fips, counties_shapefile)
            if county_geom !== nothing
                idx = GetIndex(county_geom, grid)
                locIndex[key] = idx
            else
                locIndex[key] = IndexInfo(Int[], Int[], Float64[], false, false)
            end
        else
            locIndex[key] = IndexInfo(Int[], Int[], Float64[], false, false)
        end
    end

    return locIndex
end

"""
    refine_indices_with_surrogates(locIndex::Dict{String, IndexInfo},
        county_surrogates::Dict{String, SparseMatrixCSC{Float64, Int}}) -> Dict{String, IndexInfo}

Refine area-source location indices using pre-computed county surrogate matrices.

For each FIPS-keyed entry in `locIndex` that has a matching entry in
`county_surrogates`, replaces the area-based grid allocation with
surrogate-weighted fractions via [`update_locIndex`](@ref).

Point source entries (keys containing coordinates) are left unchanged,
since point sources are allocated to their containing grid cell directly.

# Arguments
- `locIndex`: Location index dictionary from [`compute_grid_indices`](@ref).
- `county_surrogates`: Pre-computed surrogate matrices from
  [`generate_countySurrogate`](@ref), mapping FIPS codes to normalized
  allocation matrices.
"""
function refine_indices_with_surrogates(
        locIndex::Dict{String, IndexInfo},
        county_surrogates::Dict{String, SparseMatrixCSC{Float64, Int}}
    )
    for key in collect(keys(locIndex))
        # Only refine area sources (FIPS-only keys without coordinate suffix)
        if !contains(key, '_') && haskey(county_surrogates, key)
            update_locIndex(locIndex, key, county_surrogates)
        end
    end
    return locIndex
end

"""
    allocate_emissions_to_grid(emissions::DataFrame,
        locIndex::Dict{String, IndexInfo}, grid::GridDef;
        pollutant_groups::Vector{String}=["VOC", "NOX", "NH3", "SO2", "PM25"]) -> DataFrame

Distribute emissions records to grid cells using pre-computed location indices.

For each emissions record:
1. Looks up the grid cell allocation from `locIndex` using [`location_key`](@ref)
2. Distributes `ANN_VALUE` across grid cells according to fractional coverage
3. Accumulates results by `(cellIndex, SCC)` combination

# Arguments
- `emissions::DataFrame`: Must have columns `:FIPS`, `:SCC`, `:POLID`, `:ANN_VALUE`.
  May have `:LONGITUDE` and `:LATITUDE` for point sources. `ANN_VALUE` may have
  Unitful units (they will be stripped).
- `locIndex`: Location index from [`compute_grid_indices`](@ref) or
  [`refine_indices_with_surrogates`](@ref).
- `grid::GridDef`: Grid definition (used for cell index computation).
- `pollutant_groups`: Pollutant names to include in output columns.

# Returns
A `DataFrame` with columns: `cellIndex`, one column per pollutant group, and `SCC`.
Each row represents accumulated emissions for a unique `(cellIndex, SCC)` combination,
sorted by `cellIndex`.
"""
function allocate_emissions_to_grid(
        emissions::DataFrame,
        locIndex::Dict{String, IndexInfo}, grid::GridDef;
        pollutant_groups::Vector{String} = ["VOC", "NOX", "NH3", "SO2", "PM25"]
    )
    has_lon = hasproperty(emissions, :LONGITUDE)
    has_lat = hasproperty(emissions, :LATITUDE)

    # Accumulator: (cell_index, SCC) => Dict(pollutant => emission_sum)
    accumulator = Dict{Tuple{Int, String}, Dict{String, Float64}}()

    for row in eachrow(emissions)
        polid = string(row.POLID)
        polid in pollutant_groups || continue

        fips = string(row.FIPS)
        lon = has_lon ? row.LONGITUDE : missing
        lat = has_lat ? row.LATITUDE : missing
        key = location_key(fips, lon, lat)

        haskey(locIndex, key) || continue
        idx = locIndex[key]
        idx.inGrid || continue

        ann_value = Float64(ustrip(row.ANN_VALUE))
        scc = string(row.SCC)

        for k in eachindex(idx.rows)
            cell_idx = (idx.rows[k] - 1) * grid.Nx + idx.cols[k]
            acc_key = (cell_idx, scc)

            if !haskey(accumulator, acc_key)
                accumulator[acc_key] = Dict(p => 0.0 for p in pollutant_groups)
            end

            accumulator[acc_key][polid] += ann_value * idx.fracs[k]
        end
    end

    # Convert accumulator to DataFrame
    n = length(accumulator)
    cell_indices = Vector{Int}(undef, n)
    sccs = Vector{String}(undef, n)
    poll_cols = Dict(p => Vector{Float64}(undef, n) for p in pollutant_groups)

    for (i, ((cell_idx, scc), poll_values)) in enumerate(accumulator)
        cell_indices[i] = cell_idx
        sccs[i] = scc
        for p in pollutant_groups
            poll_cols[p][i] = poll_values[p]
        end
    end

    result = DataFrame(:cellIndex => cell_indices, :SCC => sccs)
    for p in pollutant_groups
        result[!, Symbol(p)] = poll_cols[p]
    end

    sort!(result, :cellIndex)
    return result
end

"""
    process_emissions_spatial(emissions::DataFrame, grid::GridDef;
        counties_shapefile::String="",
        county_surrogates::Dict{String, SparseMatrixCSC{Float64, Int}}=Dict{String, SparseMatrixCSC{Float64, Int}}(),
        pollutant_groups::Vector{String}=["VOC", "NOX", "NH3", "SO2", "PM25"]) -> DataFrame

Execute the complete spatial processing workflow on prepared emissions data.

This function orchestrates the full spatial allocation pipeline:
1. Compute grid cell indices for all unique locations ([`compute_grid_indices`](@ref))
2. Optionally refine area-source indices with surrogate data ([`refine_indices_with_surrogates`](@ref))
3. Allocate emissions to grid cells ([`allocate_emissions_to_grid`](@ref))

The input `emissions` DataFrame should already be processed through the data
preparation pipeline ([`filter_known_pollutants`](@ref), [`map_pollutant_names!`](@ref)).

# Arguments
- `emissions::DataFrame`: Prepared emissions data with columns `:FIPS`, `:SCC`,
  `:POLID`, `:ANN_VALUE`. May have `:LONGITUDE`/`:LATITUDE` for point sources.
- `grid::GridDef`: Target grid definition.
- `counties_shapefile::String=""`: Path to counties shapefile for area source indexing.
- `county_surrogates`: Pre-computed surrogate matrices for index refinement.
- `pollutant_groups`: Pollutant groups for output columns.

# Returns
A `DataFrame` with gridded emissions, one row per `(cellIndex, SCC)` combination.

# Example
```julia
grid = NewGridRegular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
emissions = DataFrame(
    FIPS=["00001", "00001"], POLID=["NOX", "VOC"],
    SCC=["2103007000", "2103007000"],
    ANN_VALUE=[1.0e-3, 5.0e-4],
    LONGITUDE=[0.5, 0.5], LATITUDE=[0.5, 0.5])
result = process_emissions_spatial(emissions, grid)
```
"""
function process_emissions_spatial(
        emissions::DataFrame, grid::GridDef;
        counties_shapefile::String = "",
        county_surrogates::Dict{String, SparseMatrixCSC{Float64, Int}} = Dict{String, SparseMatrixCSC{Float64, Int}}(),
        pollutant_groups::Vector{String} = ["VOC", "NOX", "NH3", "SO2", "PM25"]
    )
    # Step 1: Compute grid indices for all unique locations
    locIndex = compute_grid_indices(emissions, grid; counties_shapefile = counties_shapefile)

    # Step 2: Refine with surrogates if provided
    if !isempty(county_surrogates)
        refine_indices_with_surrogates(locIndex, county_surrogates)
    end

    # Step 3: Allocate emissions to grid cells
    return allocate_emissions_to_grid(
        emissions, locIndex, grid;
        pollutant_groups = pollutant_groups
    )
end
