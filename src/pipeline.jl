export read_ff10, aggregate_emissions, filter_known_pollutants,
    map_pollutant_names!, assign_surrogates, build_data_weight_map,
    process_emissions

"""
    read_ff10(filepath::AbstractString, format::Symbol) -> EmissionsDataFrame

Read an FF10-format CSV file and return the appropriate `EmissionsDataFrame` subtype.

`format` must be one of `:nonpoint`, `:point`, `:nonroad`, or `:onroad`.

# Examples
```julia
emis = read_ff10("nonpoint_2019.csv", :nonpoint)
emis.df  # access the underlying DataFrame
```
"""
function read_ff10(filepath::AbstractString, format::Symbol)
    df = CSV.read(filepath, DataFrame; header = false, comment = "#", silencewarnings = true)
    if format == :nonpoint
        return FF10NonPointDataFrame(df)
    elseif format == :point
        return FF10PointDataFrame(df)
    elseif format == :nonroad
        return FF10NonRoadDataFrame(df)
    elseif format == :onroad
        return FF10OnRoadDataFrame(df)
    else
        throw(ArgumentError("Unknown FF10 format: $format. Must be :nonpoint, :point, :nonroad, or :onroad."))
    end
end

"""
    aggregate_emissions(dfs::Vector{DataFrame}) -> DataFrame

Concatenate multiple emission DataFrames, adding `LONGITUDE` and `LATITUDE` columns
(as `missing`) for area sources that lack them, then group by
`[:POLID, :COUNTRY, :FIPS, :SCC, :LONGITUDE, :LATITUDE]` and sum `ANN_VALUE`.

The result is sorted by `ANN_VALUE` descending.
"""
function aggregate_emissions(dfs::Vector{DataFrame})
    combined_dfs = DataFrame[]
    for df in dfs
        d = copy(df)
        if !hasproperty(d, :LONGITUDE)
            d[!, :LONGITUDE] = fill(missing, nrow(d))
        end
        if !hasproperty(d, :LATITUDE)
            d[!, :LATITUDE] = fill(missing, nrow(d))
        end
        push!(combined_dfs, select(d, :POLID, :COUNTRY, :FIPS, :SCC, :LONGITUDE, :LATITUDE, :ANN_VALUE))
    end
    combined = vcat(combined_dfs...)
    grouped = combine(
        groupby(combined, [:POLID, :COUNTRY, :FIPS, :SCC, :LONGITUDE, :LATITUDE]),
        :ANN_VALUE => sum => :ANN_VALUE
    )
    sort!(grouped, :ANN_VALUE, rev = true)
    return grouped
end

"""
    filter_known_pollutants(df::DataFrame) -> DataFrame

Keep only rows where `POLID` is a key in the [`Pollutants`](@ref) dictionary.
Returns a filtered copy.
"""
function filter_known_pollutants(df::DataFrame)
    return filter(row -> haskey(Pollutants, row.POLID), df)
end

"""
    map_pollutant_names!(df::DataFrame) -> DataFrame

Map `POLID` values to standard pollutant group names using the [`Pollutants`](@ref) dictionary (in-place).
For example, "EXH__VOC" becomes "VOC" and "PM25-PRI" becomes "PM25".
"""
function map_pollutant_names!(df::DataFrame)
    df[!, :POLID] = [get(Pollutants, p, p) for p in df[!, :POLID]]
    return df
end

"""
    assign_surrogates(emissions::DataFrame, gridref::DataFrame) -> DataFrame

Assign spatial surrogates to emissions records by joining with grid reference data.

Performs a left join on `[:COUNTRY, :FIPS, :SCC]`. For rows that remain unmatched
(no surrogate assigned), falls back to a state-level match by retrying the join
with `FIPS="00000"` in the grid reference, then restoring the original FIPS code.

Returns a copy of `emissions` with a `Surrogate` column added.
"""
function assign_surrogates(emissions::DataFrame, gridref::DataFrame)
    emis = copy(emissions)

    # First pass: direct match on COUNTRY + FIPS + SCC
    joined = leftjoin(emis, gridref; on = [:COUNTRY, :FIPS, :SCC], matchmissing = :notequal)

    # Identify unmatched rows
    unmatched_mask = ismissing.(joined.Surrogate)
    if any(unmatched_mask)
        # Build fallback gridref with FIPS="00000" entries
        fallback_ref = filter(row -> row.FIPS == "00000", gridref)
        if nrow(fallback_ref) > 0
            fallback_ref = select(fallback_ref, :COUNTRY, :SCC, :Surrogate)

            unmatched = joined[unmatched_mask, :]
            select!(unmatched, Not(:Surrogate))
            original_fips = unmatched[!, :FIPS]

            # Join unmatched on COUNTRY + SCC only (state-level fallback)
            fallback_joined = leftjoin(unmatched, fallback_ref; on = [:COUNTRY, :SCC], matchmissing = :notequal)
            # Restore original FIPS
            fallback_joined[!, :FIPS] = original_fips

            # Replace unmatched rows in the result
            joined[unmatched_mask, :] = fallback_joined
        end
    end

    return joined
end

"""
    build_data_weight_map(emissions::DataFrame, srgSpecs::Vector{SurrogateSpec}) -> Dict{Tuple{String,String}, Vector{String}}

Identify which surrogates are actually needed based on the emissions data, and build a
mapping from `(data_shapefile, weight_shapefile)` pairs to lists of `"Region+Code"` strings.

This avoids loading shapefiles that aren't referenced by any emissions record.
"""
function build_data_weight_map(emissions::DataFrame, srgSpecs::Vector{SurrogateSpec})
    result = Dict{Tuple{String, String}, Vector{String}}()

    # Get unique (COUNTRY, Surrogate) pairs from emissions
    if !hasproperty(emissions, :Surrogate)
        return result
    end

    for row in eachrow(emissions)
        ismissing(row.Surrogate) && continue
        region = row.COUNTRY
        code = row.Surrogate
        srg = find_surrogate_by_code(srgSpecs, region, code)
        if srg === nothing
            # Try without region filter
            srg = find_surrogate_by_code(srgSpecs, code)
        end
        srg === nothing && continue

        key = (srg.DataShapefile, srg.WeightShapefile)
        label = "$(srg.Region)+$(srg.Code)"
        if !haskey(result, key)
            result[key] = String[]
        end
        if !(label in result[key])
            push!(result[key], label)
        end
    end

    return result
end

"""
    process_emissions(;
        inventory_files::Vector{Tuple{String, Symbol}},
        grid::GridDef,
        profiles::DataFrame=DataFrame(),
        xref::DataFrame=DataFrame(),
        episode_start::DateTime=DateTime(2019, 1, 1),
        episode_end::DateTime=DateTime(2019, 1, 2),
        gridref::DataFrame=DataFrame(),
        county_surrogates::Dict{String, SparseMatrixCSC{Float64, Int}}=Dict{String, SparseMatrixCSC{Float64, Int}}(),
        counties_shapefile::String="",
        pollutant_groups::Vector{String}=["VOC", "NOX", "NH3", "SO2", "PM25"],
        do_temporal::Bool=true,
        timezone_offset::Int=0
    ) -> DataFrame

Execute the complete emissions processing pipeline (all in-memory).

Pipeline steps:
1. Read FF10 inventory files
2. Aggregate and filter emissions
3. Map pollutant names
4. Optionally assign surrogates (if gridref provided)
5. Spatial allocation to grid
6. Temporal allocation (if profiles and xref provided and `do_temporal=true`)
7. Merge spatial + temporal into gridded hourly emissions

# Arguments
- `inventory_files`: Vector of `(filepath, format)` tuples for FF10 files.
- `grid`: Target grid definition.
- `profiles`: Temporal profiles from [`read_temporal_profiles`](@ref).
- `xref`: Temporal cross-reference from [`read_temporal_xref`](@ref).
- `episode_start`, `episode_end`: Episode time range for temporal allocation.
- `gridref`: Grid reference DataFrame for surrogate assignment.
- `county_surrogates`: Pre-computed surrogate matrices.
- `counties_shapefile`: Path to counties shapefile.
- `pollutant_groups`: Pollutant groups to process.
- `do_temporal`: Whether to perform temporal allocation.
- `timezone_offset`: UTC offset in hours.

# Returns
If `do_temporal=true`: A `DataFrame` with gridded hourly emissions
(columns: `:grid_row`, `:grid_col`, `:hour`, `:pollutant`, `:emission_rate`).
If `do_temporal=false`: A `DataFrame` with annual gridded emissions from
[`process_emissions_spatial`](@ref).
"""
function process_emissions(;
        inventory_files::Vector{Tuple{String, Symbol}},
        grid::GridDef,
        profiles::DataFrame = DataFrame(),
        xref::DataFrame = DataFrame(),
        episode_start::DateTime = DateTime(2019, 1, 1),
        episode_end::DateTime = DateTime(2019, 1, 2),
        gridref::DataFrame = DataFrame(),
        county_surrogates::Dict{String, SparseMatrixCSC{Float64, Int}} = Dict{String, SparseMatrixCSC{Float64, Int}}(),
        counties_shapefile::String = "",
        pollutant_groups::Vector{String} = ["VOC", "NOX", "NH3", "SO2", "PM25"],
        do_temporal::Bool = true,
        timezone_offset::Int = 0
    )
    # Step 1: Read inventory files
    raw_dfs = DataFrame[]
    for (filepath, format) in inventory_files
        emis = read_ff10(filepath, format)
        push!(raw_dfs, emis.df)
    end

    # Step 2: Aggregate emissions
    emissions = aggregate_emissions(raw_dfs)

    # Step 3: Filter and map pollutant names
    emissions = filter_known_pollutants(emissions)
    map_pollutant_names!(emissions)

    # Step 4: Assign surrogates if gridref provided
    if nrow(gridref) > 0
        emissions = assign_surrogates(emissions, gridref)
    end

    # Step 5: Spatial allocation
    gridded = process_emissions_spatial(
        emissions, grid;
        counties_shapefile = counties_shapefile,
        county_surrogates = county_surrogates,
        pollutant_groups = pollutant_groups
    )

    if !do_temporal
        return gridded
    end

    # Step 6: Temporal allocation
    if nrow(profiles) > 0 && nrow(xref) > 0
        hourly = temporal_allocate(
            emissions, profiles, xref,
            episode_start, episode_end; timezone_offset = timezone_offset
        )

        # Step 7: Merge spatial + temporal
        locIndex = compute_grid_indices(
            emissions, grid;
            counties_shapefile = counties_shapefile
        )
        if !isempty(county_surrogates)
            refine_indices_with_surrogates(locIndex, county_surrogates)
        end

        return merge_emissions(
            hourly, locIndex, grid;
            pollutant_groups = pollutant_groups
        )
    end

    return gridded
end
