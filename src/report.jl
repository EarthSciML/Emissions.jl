export ReportConfig, summary_by_pollutant, summary_by_region, summary_by_scc,
    summary_by_time, compare_inventories, emissions_report

"""
    ReportConfig

Configuration for emissions summary reports.

# Fields
- `group_by::Vector{Symbol}`: Columns to group by (e.g., `[:FIPS, :POLID]`).
- `time_resolution::Symbol`: Temporal aggregation level (`:hourly`, `:daily`, `:monthly`, `:annual`).
- `top_n::Int`: Number of top entries to include (0 for all).
"""
struct ReportConfig
    group_by::Vector{Symbol}
    time_resolution::Symbol
    top_n::Int

    function ReportConfig(;
            group_by::Vector{Symbol} = [:POLID],
            time_resolution::Symbol = :annual,
            top_n::Int = 0)
        valid_resolutions = (:hourly, :daily, :monthly, :annual)
        if time_resolution ∉ valid_resolutions
            throw(ArgumentError(
                "time_resolution must be one of $valid_resolutions, got :$time_resolution"
            ))
        end
        return new(group_by, time_resolution, top_n)
    end
end

"""
    summary_by_pollutant(emissions::DataFrame) -> DataFrame

Compute total, mean, and maximum emission values per pollutant.

Expects a DataFrame with at least `:POLID` and `:ANN_VALUE` columns (pre-merge inventory)
or `:pollutant` and `:emission_rate` columns (post-merge gridded data).

Returns a DataFrame with columns: `:pollutant`, `:total`, `:mean`, `:max`, `:count`.
"""
function summary_by_pollutant(emissions::DataFrame)
    if nrow(emissions) == 0
        return DataFrame(
            pollutant = String[], total = Float64[], mean = Float64[],
            max = Float64[], count = Int[]
        )
    end

    # Detect column names (pre-merge vs post-merge)
    if hasproperty(emissions, :pollutant) && hasproperty(emissions, :emission_rate)
        pol_col = :pollutant
        val_col = :emission_rate
    elseif hasproperty(emissions, :POLID) && hasproperty(emissions, :ANN_VALUE)
        pol_col = :POLID
        val_col = :ANN_VALUE
    else
        throw(ArgumentError(
            "DataFrame must have (:pollutant, :emission_rate) or (:POLID, :ANN_VALUE) columns"
        ))
    end

    gdf = groupby(emissions, pol_col)
    result = combine(gdf,
        val_col => (v -> sum(Float64.(ustrip.(v)))) => :total,
        val_col => (v -> sum(Float64.(ustrip.(v))) / length(v)) => :mean,
        val_col => (v -> maximum(Float64.(ustrip.(v)))) => :max,
        nrow => :count
    )
    rename!(result, pol_col => :pollutant)
    sort!(result, :total, rev = true)
    return result
end

"""
    summary_by_region(emissions::DataFrame,
        region_map::Dict{Tuple{Int,Int}, String}) -> DataFrame

Compute total emissions per region and pollutant for gridded data.

`region_map` maps `(grid_row, grid_col)` to a region name string.

Returns a DataFrame with columns: `:region`, `:pollutant`, `:total`.
"""
function summary_by_region(emissions::DataFrame,
        region_map::Dict{Tuple{Int, Int}, String})
    if nrow(emissions) == 0
        return DataFrame(region = String[], pollutant = String[], total = Float64[])
    end

    # Add region column
    df = copy(emissions)
    df[!, :region] = [get(region_map, (r.grid_row, r.grid_col), "Unknown")
                      for r in eachrow(df)]

    gdf = groupby(df, [:region, :pollutant])
    result = combine(gdf,
        :emission_rate => sum => :total
    )
    sort!(result, :total, rev = true)
    return result
end

"""
    summary_by_scc(emissions::DataFrame) -> DataFrame

Compute total emissions per SCC code and pollutant from pre-merge inventory data.

Expects columns `:SCC`, `:POLID`, `:ANN_VALUE`.

Returns a DataFrame with columns: `:SCC`, `:POLID`, `:total`, `:count`,
sorted by total descending.
"""
function summary_by_scc(emissions::DataFrame)
    if nrow(emissions) == 0
        return DataFrame(
            SCC = String[], POLID = String[], total = Float64[], count = Int[]
        )
    end

    gdf = groupby(emissions, [:SCC, :POLID])
    result = combine(gdf,
        :ANN_VALUE => (v -> sum(Float64.(ustrip.(v)))) => :total,
        nrow => :count
    )
    sort!(result, :total, rev = true)
    return result
end

"""
    summary_by_time(emissions::DataFrame; resolution::Symbol=:daily) -> DataFrame

Aggregate gridded hourly emissions by time period.

`resolution` must be `:hourly`, `:daily`, or `:monthly`.

Returns a DataFrame with columns: `:period`, `:pollutant`, `:total`.
"""
function summary_by_time(emissions::DataFrame; resolution::Symbol = :daily)
    if nrow(emissions) == 0
        return DataFrame(period = String[], pollutant = String[], total = Float64[])
    end

    df = copy(emissions)

    if resolution == :hourly
        df[!, :period] = string.(df.hour)
    elseif resolution == :daily
        df[!, :period] = [string(Date(h)) for h in df.hour]
    elseif resolution == :monthly
        df[!, :period] = [Dates.format(h, "yyyy-mm") for h in df.hour]
    else
        throw(ArgumentError("resolution must be :hourly, :daily, or :monthly"))
    end

    gdf = groupby(df, [:period, :pollutant])
    result = combine(gdf, :emission_rate => sum => :total)
    sort!(result, [:period, :pollutant])
    return result
end

"""
    compare_inventories(base::DataFrame, scenario::DataFrame;
        join_on::Vector{Symbol}=[:FIPS, :SCC, :POLID]) -> DataFrame

Compare two emission inventories and compute absolute and percent differences.

Both DataFrames must have `join_on` columns plus `:ANN_VALUE`.

Returns a DataFrame with columns from `join_on` plus:
`:base_value`, `:scenario_value`, `:abs_diff`, `:pct_diff`.
"""
function compare_inventories(base::DataFrame, scenario::DataFrame;
        join_on::Vector{Symbol} = [:FIPS, :SCC, :POLID])
    if nrow(base) == 0 && nrow(scenario) == 0
        cols = vcat(join_on, [:base_value, :scenario_value, :abs_diff, :pct_diff])
        types = vcat(fill(String, length(join_on)), fill(Float64, 4))
        return DataFrame([c => t[] for (c, t) in zip(cols, types)])
    end

    b = select(base, join_on..., :ANN_VALUE)
    s = select(scenario, join_on..., :ANN_VALUE)

    # Aggregate each by join keys
    b_agg = combine(groupby(b, join_on),
        :ANN_VALUE => (v -> sum(Float64.(ustrip.(v)))) => :base_value)
    s_agg = combine(groupby(s, join_on),
        :ANN_VALUE => (v -> sum(Float64.(ustrip.(v)))) => :scenario_value)

    result = outerjoin(b_agg, s_agg; on = join_on, matchmissing = :equal)

    # Fill missing with 0
    result.base_value = coalesce.(result.base_value, 0.0)
    result.scenario_value = coalesce.(result.scenario_value, 0.0)

    result[!, :abs_diff] = result.scenario_value .- result.base_value
    result[!, :pct_diff] = [
        bv == 0.0 ? (sv == 0.0 ? 0.0 : Inf) : 100.0 * (sv - bv) / abs(bv)
        for (bv, sv) in zip(result.base_value, result.scenario_value)
    ]

    sort!(result, :abs_diff, rev = true, by = abs)
    return result
end

"""
    emissions_report(emissions::DataFrame; config::ReportConfig=ReportConfig()) -> DataFrame

Generate a summary report grouped by the specified columns.

Uses `config.group_by` for grouping, applies `config.top_n` limit if > 0.
Works with both pre-merge (`:POLID`, `:ANN_VALUE`) and post-merge
(`:pollutant`, `:emission_rate`) DataFrames.
"""
function emissions_report(emissions::DataFrame;
        config::ReportConfig = ReportConfig())
    if nrow(emissions) == 0
        return DataFrame()
    end

    # Detect value column
    val_col = hasproperty(emissions, :emission_rate) ? :emission_rate : :ANN_VALUE

    # Ensure all group_by columns exist
    valid_cols = filter(c -> hasproperty(emissions, c), config.group_by)
    if isempty(valid_cols)
        return DataFrame()
    end

    gdf = groupby(emissions, valid_cols)
    result = combine(gdf,
        val_col => (v -> sum(Float64.(ustrip.(v)))) => :total,
        val_col => (v -> sum(Float64.(ustrip.(v))) / length(v)) => :mean,
        nrow => :count
    )
    sort!(result, :total, rev = true)

    if config.top_n > 0 && nrow(result) > config.top_n
        result = first(result, config.top_n)
    end

    return result
end
