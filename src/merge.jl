export merge_emissions, merge_categories

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
