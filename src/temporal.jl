export read_temporal_profiles, read_temporal_xref, temporal_allocate

"""
    read_temporal_profiles(filepath::AbstractString) -> DataFrame

Read a SMOKE-format temporal profile file (ATPRO/MTPRO).

The file contains monthly, weekly, and diurnal profile factors.
Lines starting with `#` or `/PROFILE/` are skipped.
Each data line has the format: `profile_type, profile_id, factor1, factor2, ...`

Profile types:
- `MONTHLY`: 12 factors (Jan-Dec) that sum to 1.0
- `WEEKLY`: 7 factors (Mon-Sun) that sum to 7.0
- `DIURNAL` or `ALLDAY`: 24 factors (hours 0-23) that sum to 1.0
- Day-specific: `MONDAY`, `TUESDAY`, ..., `SUNDAY`, `WEEKDAY`, `WEEKEND`

Returns a DataFrame with columns: `[:profile_type, :profile_id, :factors]`
where `factors` is a `Vector{Float64}`.
"""
function read_temporal_profiles(filepath::AbstractString)
    lines = readlines(filepath)
    records = DataFrame(
        profile_type = String[],
        profile_id = Int[],
        factors = Vector{Float64}[]
    )
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#") || startswith(line, "/")
            continue
        end
        parts = split(line, r"[,\s]+")
        if length(parts) < 3
            continue
        end
        profile_type = uppercase(strip(parts[1]))
        profile_id = parse(Int, strip(parts[2]))
        factors = Float64[]
        for i in 3:length(parts)
            s = strip(parts[i])
            if !isempty(s)
                try
                    push!(factors, parse(Float64, s))
                catch
                    break
                end
            end
        end
        if !isempty(factors)
            push!(records, (profile_type = profile_type, profile_id = profile_id, factors = factors))
        end
    end
    return records
end

"""
    read_temporal_xref(filepath::AbstractString) -> DataFrame

Read a SMOKE-format temporal cross-reference file (ATREF/MTREF).

Maps SCC codes to temporal profile IDs. Lines starting with `#` are skipped.
Each data line has the format:
`FIPS;SCC;monthly_profile_id;weekly_profile_id;diurnal_profile_id[!comment]`

Returns a DataFrame with columns:
`[:FIPS, :SCC, :monthly_id, :weekly_id, :diurnal_id]`
"""
function read_temporal_xref(filepath::AbstractString)
    lines = readlines(filepath)
    records = DataFrame(
        FIPS = String[],
        SCC = String[],
        monthly_id = Int[],
        weekly_id = Int[],
        diurnal_id = Int[]
    )
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        # Remove inline comments
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ';')
        if length(parts) >= 5
            fips_raw = strip(parts[1])
            scc = strip(parts[2])
            monthly_id = parse(Int, strip(parts[3]))
            weekly_id = parse(Int, strip(parts[4]))
            diurnal_id = parse(Int, strip(parts[5]))

            # Normalize FIPS
            if length(fips_raw) == 6
                fips = lpad(fips_raw[2:end], 5, '0')
            else
                fips = lpad(fips_raw, 5, '0')
            end

            push!(
                records, (
                    FIPS = fips, SCC = scc, monthly_id = monthly_id,
                    weekly_id = weekly_id, diurnal_id = diurnal_id,
                )
            )
        end
    end
    return records
end

# Pre-built Dict cache type for profile lookup
const _ProfileKey = Tuple{String, Int}
const _ProfileCache = Dict{_ProfileKey, Vector{Float64}}

"""
    _build_profile_cache(profiles::DataFrame) -> Dict

Build a Dict-based lookup cache from profiles DataFrame for O(1) lookup.
"""
function _build_profile_cache(profiles::DataFrame)
    cache = _ProfileCache()
    for row in eachrow(profiles)
        cache[(row.profile_type, row.profile_id)] = row.factors
    end
    return cache
end

"""
    _lookup_profile(profiles::DataFrame, profile_type::String, profile_id::Int) -> Vector{Float64}

Internal helper to look up a temporal profile by type and ID.
Returns the factor vector or a uniform default if not found.
"""
function _lookup_profile(profiles::DataFrame, profile_type::String, profile_id::Int)
    for row in eachrow(profiles)
        if row.profile_type == profile_type && row.profile_id == profile_id
            return row.factors
        end
    end
    return _default_profile(profile_type)
end

"""
    _lookup_profile_cached(cache::Dict, profile_type::String, profile_id::Int) -> Vector{Float64}

O(1) profile lookup using pre-built cache.
"""
function _lookup_profile_cached(cache::_ProfileCache, profile_type::String, profile_id::Int)
    return get(cache, (profile_type, profile_id), _default_profile(profile_type))
end

"""
    _default_profile(profile_type::String) -> Vector{Float64}

Return uniform default profile factors for the given profile type.
"""
function _default_profile(profile_type::String)
    if profile_type == "MONTHLY"
        return fill(1.0 / 12.0, 12)
    elseif profile_type == "WEEKLY"
        return fill(1.0, 7)
    else  # DIURNAL / ALLDAY / day-specific
        return fill(1.0 / 24.0, 24)
    end
end

# Day name constants for day-specific diurnal profile matching
const _DOW_NAMES = ["MONDAY", "TUESDAY", "WEDNESDAY", "THURSDAY", "FRIDAY", "SATURDAY", "SUNDAY"]

"""
    _lookup_diurnal_profile(cache::Dict, profile_id::Int, dow::Int) -> Vector{Float64}

Look up a diurnal profile with SMOKE-style day-specific fallback hierarchy:
1. DIURNAL (exact profile)
2. Day-specific name (MONDAY, TUESDAY, etc.)
3. WEEKDAY (Mon-Fri) or WEEKEND (Sat-Sun)
4. ALLDAY
5. Uniform default
"""
function _lookup_diurnal_profile(cache::_ProfileCache, profile_id::Int, dow::Int)
    # Try DIURNAL first
    result = get(cache, ("DIURNAL", profile_id), nothing)
    if result !== nothing && length(result) == 24
        return result
    end

    # Try day-specific name
    day_name = _DOW_NAMES[dow]
    result = get(cache, (day_name, profile_id), nothing)
    if result !== nothing && length(result) == 24
        return result
    end

    # Try WEEKDAY/WEEKEND
    if dow <= 5
        result = get(cache, ("WEEKDAY", profile_id), nothing)
    else
        result = get(cache, ("WEEKEND", profile_id), nothing)
    end
    if result !== nothing && length(result) == 24
        return result
    end

    # Try ALLDAY
    result = get(cache, ("ALLDAY", profile_id), nothing)
    if result !== nothing && length(result) == 24
        return result
    end

    return fill(1.0 / 24.0, 24)
end

# Pre-built xref cache type
const _XrefKey = Tuple{String, String}
const _XrefValue = NamedTuple{(:monthly_id, :weekly_id, :diurnal_id), Tuple{Int, Int, Int}}
const _XrefCache = Dict{_XrefKey, _XrefValue}

"""
    _build_xref_cache(xref::DataFrame) -> Dict

Build a Dict-based lookup cache from xref DataFrame for O(1) lookup.
"""
function _build_xref_cache(xref::DataFrame)
    cache = _XrefCache()
    for row in eachrow(xref)
        cache[(row.FIPS, row.SCC)] = (monthly_id = row.monthly_id, weekly_id = row.weekly_id, diurnal_id = row.diurnal_id)
    end
    return cache
end

"""
    _match_temporal_xref(xref::DataFrame, fips::AbstractString, scc::AbstractString) -> NamedTuple

Match an emissions record to temporal profile IDs using hierarchical matching:
1. Exact FIPS + SCC match
2. SCC match with FIPS="00000" (national default)
3. Default profile IDs (1, 1, 1) if no match found
"""
function _match_temporal_xref(xref::DataFrame, fips::AbstractString, scc::AbstractString)
    # Try exact FIPS + SCC match
    for row in eachrow(xref)
        if row.FIPS == fips && row.SCC == scc
            return (monthly_id = row.monthly_id, weekly_id = row.weekly_id, diurnal_id = row.diurnal_id)
        end
    end
    # Try national default (FIPS="00000") + SCC match
    for row in eachrow(xref)
        if row.FIPS == "00000" && row.SCC == scc
            return (monthly_id = row.monthly_id, weekly_id = row.weekly_id, diurnal_id = row.diurnal_id)
        end
    end
    # Try SCC-only match (any FIPS)
    for row in eachrow(xref)
        if row.SCC == scc
            return (monthly_id = row.monthly_id, weekly_id = row.weekly_id, diurnal_id = row.diurnal_id)
        end
    end
    # Default
    return (monthly_id = 1, weekly_id = 1, diurnal_id = 1)
end

"""
    _match_temporal_xref_cached(cache::Dict, fips::AbstractString, scc::AbstractString) -> NamedTuple

O(1) temporal xref matching using pre-built cache with hierarchical fallback.
"""
function _match_temporal_xref_cached(cache::_XrefCache, fips::AbstractString, scc::AbstractString)
    # Level 1: exact FIPS + SCC
    result = get(cache, (fips, scc), nothing)
    result !== nothing && return result

    # Level 2: national default + SCC
    result = get(cache, ("00000", scc), nothing)
    result !== nothing && return result

    # Default
    return (monthly_id = 1, weekly_id = 1, diurnal_id = 1)
end

"""
    temporal_allocate(emissions::DataFrame, profiles::DataFrame,
        xref::DataFrame, episode_start::DateTime, episode_end::DateTime;
        timezone_offset::Int=0,
        timezone_map::Dict{String,Int}=Dict{String,Int}()) -> DataFrame

Convert annual emissions to hourly emissions using temporal profiles.

Applies monthly, day-of-week, and diurnal temporal factors to distribute
annual emissions to hourly values for the specified episode period.

# Arguments
- `emissions::DataFrame`: Annual emissions with columns `:FIPS`, `:SCC`, `:POLID`, `:ANN_VALUE`.
  `ANN_VALUE` should be in mass/time units (e.g., kg/s). May also have
  `:LONGITUDE`, `:LATITUDE`, `:COUNTRY`, and `:Surrogate`.
- `profiles::DataFrame`: Temporal profiles from [`read_temporal_profiles`](@ref).
- `xref::DataFrame`: Temporal cross-reference from [`read_temporal_xref`](@ref).
- `episode_start::DateTime`: Start of the episode period.
- `episode_end::DateTime`: End of the episode period.
- `timezone_offset::Int=0`: Default UTC offset in hours for local time adjustment.
- `timezone_map::Dict{String,Int}=Dict{String,Int}()`: Optional per-FIPS timezone
  offsets (FIPS code => UTC offset). When present, overrides `timezone_offset`
  for sources with matching FIPS codes.

# Returns
A `DataFrame` with columns: `:FIPS`, `:SCC`, `:POLID`, `:hour` (DateTime),
`:emission_rate` (Float64), plus any location columns from the input.

The temporal allocation follows SMOKE conventions:
- Monthly factor: fraction of annual total for the given month (sums to 1.0)
- Weekly factor: relative weight for the day of week (sums to 7.0)
- Diurnal factor: fraction of daily total for the given hour (sums to 1.0)
- Hourly rate = ANN_VALUE × monthly_factor × (weekly_factor / 7) × (diurnal_factor × 24)
- Diurnal profiles are selected with SMOKE-style day-specific fallback:
  DIURNAL → day name (MONDAY, etc.) → WEEKDAY/WEEKEND → ALLDAY → uniform
"""
function temporal_allocate(
        emissions::DataFrame,
        profiles::DataFrame,
        xref::DataFrame,
        episode_start::DateTime,
        episode_end::DateTime;
        timezone_offset::Int = 0,
        timezone_map::Dict{String, Int} = Dict{String, Int}()
    )
    # Pre-compute hours in the episode
    n_hours = Int(Dates.value(episode_end - episode_start)) ÷ 3_600_000
    hours = [episode_start + Hour(h) for h in 0:(n_hours - 1)]

    if isempty(hours) || nrow(emissions) == 0
        return DataFrame(
            FIPS = String[],
            SCC = String[],
            POLID = String[],
            hour = DateTime[],
            emission_rate = Float64[]
        )
    end

    n_sources = nrow(emissions)
    n_total = n_sources * length(hours)

    # Build caches for O(1) lookups
    profile_cache = _build_profile_cache(profiles)
    xref_cache = _build_xref_cache(xref)

    # Pre-allocate output column vectors
    out_fips = Vector{String}(undef, n_total)
    out_scc = Vector{String}(undef, n_total)
    out_polid = Vector{String}(undef, n_total)
    out_hour = Vector{DateTime}(undef, n_total)
    out_rate = Vector{Float64}(undef, n_total)

    idx = 0
    for erow in eachrow(emissions)
        fips = string(erow.FIPS)
        scc = string(erow.SCC)
        polid = string(erow.POLID)
        ann_value = Float64(ustrip(erow.ANN_VALUE))

        # Per-source timezone: use timezone_map if available, else global default
        tz_offset = get(timezone_map, fips, timezone_offset)

        # Look up temporal profile IDs
        profile_ids = _match_temporal_xref_cached(xref_cache, fips, scc)

        # Get profile factors
        monthly_factors = _lookup_profile_cached(profile_cache, "MONTHLY", profile_ids.monthly_id)
        weekly_factors = _lookup_profile_cached(profile_cache, "WEEKLY", profile_ids.weekly_id)

        for hr in hours
            local_hr = hr + Hour(tz_offset)
            month_idx = Dates.month(local_hr)
            dow = Dates.dayofweek(local_hr)  # 1=Monday, 7=Sunday
            hour_idx = Dates.hour(local_hr) + 1  # 1-indexed

            # Day-specific diurnal profile lookup
            diurnal_factors = _lookup_diurnal_profile(profile_cache, profile_ids.diurnal_id, dow)

            mf = month_idx <= length(monthly_factors) ? monthly_factors[month_idx] : 1.0 / 12.0
            wf = dow <= length(weekly_factors) ? weekly_factors[dow] : 1.0
            df = hour_idx <= length(diurnal_factors) ? diurnal_factors[hour_idx] : 1.0 / 24.0

            # hourly_rate = annual_rate * monthly_frac * (weekly_weight/7) * (diurnal_frac * 24)
            hourly_rate = ann_value * mf * (wf / 7.0) * (df * 24.0)

            idx += 1
            out_fips[idx] = fips
            out_scc[idx] = scc
            out_polid[idx] = polid
            out_hour[idx] = hr
            out_rate[idx] = hourly_rate
        end
    end

    return DataFrame(
        FIPS = out_fips,
        SCC = out_scc,
        POLID = out_polid,
        hour = out_hour,
        emission_rate = out_rate,
    )
end
