export merge_emissions, merge_categories, merge_categories_tracked,
    merge_2d_3d, to_model_ready

"""
    merge_emissions(hourly_emissions::DataFrame, locIndex::Dict{String, IndexInfo},
        grid::GridDef; pollutant_groups::Vector{String}=["VOC", "NOX", "NH3", "SO2", "PM25"]) -> DataFrame

Combine temporally-allocated hourly emissions with spatial allocation to produce
gridded hourly emissions.

For each hour and source:
1. Looks up the spatial allocation (grid cells and fractions) from `locIndex`
2. Distributes the hourly emission rate across grid cells according to fractions
3. Accumulates results by `(grid_row, grid_col, hour, pollutant)` combination

# Arguments
- `hourly_emissions::DataFrame`: From [`temporal_allocate`](@ref), with columns
  `:FIPS`, `:SCC`, `:POLID`, `:hour`, `:emission_rate`.
  May also have `:LONGITUDE` and `:LATITUDE` for point sources.
- `locIndex`: Pre-computed location index from [`compute_grid_indices`](@ref).
- `grid::GridDef`: Target grid definition.
- `pollutant_groups`: Pollutant groups to include in output.

# Returns
A `DataFrame` with columns: `:grid_row`, `:grid_col`, `:hour`, `:pollutant`, `:emission_rate`.
Each row represents accumulated emissions for a unique grid cell, hour, and pollutant.
"""
function merge_emissions(
        hourly_emissions::DataFrame,
        locIndex::Dict{String, IndexInfo},
        grid::GridDef;
        pollutant_groups::Vector{String} = ["VOC", "NOX", "NH3", "SO2", "PM25"],
        species_list::Vector{String} = String[]
    )
    # Use species_list if provided, otherwise use pollutant_groups
    filter_list = isempty(species_list) ? pollutant_groups : species_list
    has_lon = hasproperty(hourly_emissions, :LONGITUDE)
    has_lat = hasproperty(hourly_emissions, :LATITUDE)
    has_layer = hasproperty(hourly_emissions, :layer)

    if has_layer
        # Accumulator with layer: (row, col, hour, pollutant, layer) => emission_sum
        accumulator = Dict{Tuple{Int, Int, DateTime, String, Int}, Float64}()

        for row in eachrow(hourly_emissions)
            polid = string(row.POLID)
            polid in filter_list || continue

            fips = string(row.FIPS)
            lon = has_lon ? row.LONGITUDE : missing
            lat = has_lat ? row.LATITUDE : missing
            key = location_key(fips, lon, lat)

            haskey(locIndex, key) || continue
            idx = locIndex[key]
            idx.inGrid || continue

            hr = row.hour
            rate = Float64(row.emission_rate)
            layer = Int(row.layer)
            layer_frac = hasproperty(hourly_emissions, :layer_fraction) ? Float64(row.layer_fraction) : 1.0

            for k in eachindex(idx.rows)
                acc_key = (idx.rows[k], idx.cols[k], hr, polid, layer)
                accumulator[acc_key] = get(accumulator, acc_key, 0.0) + rate * idx.fracs[k] * layer_frac
            end
        end

        n = length(accumulator)
        grid_rows = Vector{Int}(undef, n)
        grid_cols = Vector{Int}(undef, n)
        hours = Vector{DateTime}(undef, n)
        pollutants = Vector{String}(undef, n)
        layers = Vector{Int}(undef, n)
        rates = Vector{Float64}(undef, n)

        for (i, ((gr, gc, hr, pol, lay), rate)) in enumerate(accumulator)
            grid_rows[i] = gr
            grid_cols[i] = gc
            hours[i] = hr
            pollutants[i] = pol
            layers[i] = lay
            rates[i] = rate
        end

        result = DataFrame(
            grid_row = grid_rows,
            grid_col = grid_cols,
            hour = hours,
            pollutant = pollutants,
            layer = layers,
            emission_rate = rates,
        )
        sort!(result, [:hour, :grid_row, :grid_col, :pollutant, :layer])
        return result
    else
        # Original accumulator without layer
        accumulator = Dict{Tuple{Int, Int, DateTime, String}, Float64}()

        for row in eachrow(hourly_emissions)
            polid = string(row.POLID)
            polid in filter_list || continue

            fips = string(row.FIPS)
            lon = has_lon ? row.LONGITUDE : missing
            lat = has_lat ? row.LATITUDE : missing
            key = location_key(fips, lon, lat)

            haskey(locIndex, key) || continue
            idx = locIndex[key]
            idx.inGrid || continue

            hr = row.hour
            rate = Float64(row.emission_rate)

            for k in eachindex(idx.rows)
                acc_key = (idx.rows[k], idx.cols[k], hr, polid)
                accumulator[acc_key] = get(accumulator, acc_key, 0.0) + rate * idx.fracs[k]
            end
        end

        n = length(accumulator)
        grid_rows = Vector{Int}(undef, n)
        grid_cols = Vector{Int}(undef, n)
        hours = Vector{DateTime}(undef, n)
        pollutants = Vector{String}(undef, n)
        rates = Vector{Float64}(undef, n)

        for (i, ((gr, gc, hr, pol), rate)) in enumerate(accumulator)
            grid_rows[i] = gr
            grid_cols[i] = gc
            hours[i] = hr
            pollutants[i] = pol
            rates[i] = rate
        end

        result = DataFrame(
            grid_row = grid_rows,
            grid_col = grid_cols,
            hour = hours,
            pollutant = pollutants,
            emission_rate = rates,
        )
        sort!(result, [:hour, :grid_row, :grid_col, :pollutant])
        return result
    end
end

"""
    merge_categories(dfs::DataFrame...) -> DataFrame
    merge_categories(dfs::Vector{DataFrame}) -> DataFrame

Combine gridded emissions from different source categories (area, point, mobile, etc.).

Stacks DataFrames and sums emissions for the same `(grid_row, grid_col, hour, pollutant)`.

# Arguments
- `dfs`: One or more DataFrames with columns `:grid_row`, `:grid_col`, `:hour`,
  `:pollutant`, `:emission_rate`.

# Returns
A combined `DataFrame` with summed emission rates, sorted by hour, row, col, pollutant.
"""
function merge_categories(dfs::DataFrame...)
    return merge_categories(collect(dfs))
end

function merge_categories(dfs::Vector{DataFrame})
    if isempty(dfs)
        return DataFrame(
            grid_row = Int[],
            grid_col = Int[],
            hour = DateTime[],
            pollutant = String[],
            emission_rate = Float64[]
        )
    end

    combined = vcat(dfs...)
    if nrow(combined) == 0
        return combined
    end

    result = combine(
        groupby(combined, [:grid_row, :grid_col, :hour, :pollutant]),
        :emission_rate => sum => :emission_rate
    )
    sort!(result, [:hour, :grid_row, :grid_col, :pollutant])
    return result
end

"""
    merge_categories_tracked(dfs::Vector{Pair{String,DataFrame}}) -> DataFrame

Combine gridded emissions from different source categories, tracking which category
each emission came from.

Each pair is `"category_name" => gridded_df`. Adds a `:source_category` column and
sums by `(grid_row, grid_col, hour, pollutant, source_category)`.

# Arguments
- `dfs`: Vector of `"name" => DataFrame` pairs.

# Returns
A combined `DataFrame` with `:source_category` column added, sorted by hour, row, col, pollutant, category.
"""
function merge_categories_tracked(dfs::Vector{Pair{String, DataFrame}})
    if isempty(dfs)
        return DataFrame(
            grid_row = Int[],
            grid_col = Int[],
            hour = DateTime[],
            pollutant = String[],
            source_category = String[],
            emission_rate = Float64[]
        )
    end

    parts = DataFrame[]
    for (name, df) in dfs
        d = copy(df)
        d[!, :source_category] .= name
        push!(parts, d)
    end

    combined = vcat(parts...)
    if nrow(combined) == 0
        return combined
    end

    result = combine(
        groupby(combined, [:grid_row, :grid_col, :hour, :pollutant, :source_category]),
        :emission_rate => sum => :emission_rate
    )
    sort!(result, [:hour, :grid_row, :grid_col, :pollutant, :source_category])
    return result
end

"""
    merge_2d_3d(surface::DataFrame, elevated::DataFrame; surface_layer::Int=1) -> DataFrame

Merge 2D surface emissions with 3D elevated emissions into a unified DataFrame
with a `:layer` column.

Surface emissions are assigned to `surface_layer` (default 1). Elevated emissions
must already have a `:layer` column.

# Returns
A `DataFrame` with columns: `:grid_row`, `:grid_col`, `:hour`, `:pollutant`,
`:layer`, `:emission_rate`.
"""
function merge_2d_3d(surface::DataFrame, elevated::DataFrame; surface_layer::Int = 1)
    # Add layer column to surface data
    surf = copy(surface)
    surf[!, :layer] .= surface_layer

    # Elevated should already have :layer
    if !hasproperty(elevated, :layer)
        elev = copy(elevated)
        elev[!, :layer] .= surface_layer
    else
        elev = elevated
    end

    # Select common columns
    cols = [:grid_row, :grid_col, :hour, :pollutant, :layer, :emission_rate]
    combined = vcat(select(surf, cols), select(elev, cols))

    if nrow(combined) == 0
        return combined
    end

    result = combine(
        groupby(combined, [:grid_row, :grid_col, :hour, :pollutant, :layer]),
        :emission_rate => sum => :emission_rate
    )
    sort!(result, [:hour, :layer, :grid_row, :grid_col, :pollutant])
    return result
end

"""
    to_model_ready(merged::DataFrame, grid::GridDef, hours::Vector{DateTime};
        n_layers::Int=1) -> Dict{String, Array{Float64}}

Convert merged gridded emissions DataFrame to dense 4D arrays suitable for
air quality model input.

# Arguments
- `merged::DataFrame`: Gridded emissions with columns `:grid_row`, `:grid_col`,
  `:hour`, `:pollutant`, `:emission_rate`, and optionally `:layer`.
- `grid::GridDef`: Grid definition (provides `Ny` rows and `Nx` columns).
- `hours::Vector{DateTime}`: Ordered list of time steps.
- `n_layers::Int=1`: Number of vertical layers.

# Returns
`Dict{String, Array{Float64, 4}}` mapping pollutant names to arrays with
dimensions `[nrows, ncols, n_layers, n_hours]`.
"""
function to_model_ready(merged::DataFrame, grid::GridDef, hours::Vector{DateTime};
        n_layers::Int = 1)
    nrows = grid.Ny
    ncols = grid.Nx
    n_hours_t = length(hours)

    # Build hour -> index mapping
    hour_idx = Dict{DateTime, Int}()
    for (i, h) in enumerate(hours)
        hour_idx[h] = i
    end

    # Get unique pollutants
    pollutants = unique(merged.pollutant)
    has_layer = hasproperty(merged, :layer)

    result = Dict{String, Array{Float64, 4}}()
    for pol in pollutants
        result[pol] = zeros(Float64, nrows, ncols, n_layers, n_hours_t)
    end

    for row in eachrow(merged)
        pol = row.pollutant
        haskey(result, pol) || continue

        r = row.grid_row
        c = row.grid_col
        h_idx = get(hour_idx, row.hour, 0)
        h_idx == 0 && continue
        l = has_layer ? row.layer : 1

        # Bounds check
        (1 <= r <= nrows && 1 <= c <= ncols && 1 <= l <= n_layers) || continue

        result[pol][r, c, l, h_idx] += row.emission_rate
    end

    return result
end
