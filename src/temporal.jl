export read_temporal_profiles, read_temporal_xref, temporal_allocate,
    read_day_specific, read_hour_specific

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
    read_day_specific(filepath::AbstractString) -> DataFrame

Read a SMOKE-format day-specific emissions file (PTDAY/ARDAY).

Each data line has the format: `FIPS;SCC;POLID;date;day_emis[!comment]`

Returns a DataFrame with columns:
`[:FIPS, :SCC, :POLID, :date (Date), :day_value (Float64)]`
"""
function read_day_specific(filepath::AbstractString)
    lines = readlines(filepath)
    records = DataFrame(
        FIPS = String[],
        SCC = String[],
        POLID = String[],
        date = Date[],
        day_value = Float64[]
    )
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ';')
        if length(parts) >= 5
            fips_raw = strip(parts[1])
            if length(fips_raw) == 6
                fips = lpad(fips_raw[2:end], 5, '0')
            else
                fips = lpad(fips_raw, 5, '0')
            end
            scc = strip(parts[2])
            polid = strip(parts[3])
            dt = Date(strip(parts[4]), dateformat"mm/dd/yyyy")
            day_val = parse(Float64, strip(parts[5]))
            push!(records, (FIPS = fips, SCC = scc, POLID = polid, date = dt, day_value = day_val))
        end
    end
    return records
end

"""
    read_hour_specific(filepath::AbstractString) -> DataFrame

Read a SMOKE-format hour-specific emissions file (PTHOURLY/ARHOURLY).

Each data line has the format: `FIPS;SCC;POLID;date;h0;h1;...;h23[!comment]`

Returns a DataFrame with columns:
`[:FIPS, :SCC, :POLID, :date (Date), :hourly_values (Vector{Float64})]`
"""
function read_hour_specific(filepath::AbstractString)
    lines = readlines(filepath)
    records = DataFrame(
        FIPS = String[],
        SCC = String[],
        POLID = String[],
        date = Date[],
        hourly_values = Vector{Float64}[]
    )
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ';')
        if length(parts) >= 28  # FIPS;SCC;POLID;date;h0..h23
            fips_raw = strip(parts[1])
            if length(fips_raw) == 6
                fips = lpad(fips_raw[2:end], 5, '0')
            else
                fips = lpad(fips_raw, 5, '0')
            end
            scc = strip(parts[2])
            polid = strip(parts[3])
            dt = Date(strip(parts[4]), dateformat"mm/dd/yyyy")
            hourly = Float64[]
            for i in 5:28
                push!(hourly, parse(Float64, strip(parts[i])))
            end
            push!(records, (FIPS = fips, SCC = scc, POLID = polid, date = dt, hourly_values = hourly))
        end
    end
    return records
end

# Cache types for day-specific and hour-specific lookups
const _DaySpecificKey = Tuple{String, String, String, Date}
const _DaySpecificCache = Dict{_DaySpecificKey, Float64}
const _HourSpecificKey = Tuple{String, String, String, Date}
const _HourSpecificCache = Dict{_HourSpecificKey, Vector{Float64}}

"""
    _build_day_specific_cache(df::DataFrame) -> Dict

Build an O(1) lookup cache from day-specific DataFrame.
Keys are `(FIPS, SCC, POLID, Date)`.
"""
function _build_day_specific_cache(df::DataFrame)
    cache = _DaySpecificCache()
    for row in eachrow(df)
        cache[(row.FIPS, row.SCC, row.POLID, row.date)] = row.day_value
    end
    return cache
end

"""
    _build_hour_specific_cache(df::DataFrame) -> Dict

Build an O(1) lookup cache from hour-specific DataFrame.
Keys are `(FIPS, SCC, POLID, Date)`.
"""
function _build_hour_specific_cache(df::DataFrame)
    cache = _HourSpecificCache()
    for row in eachrow(df)
        cache[(row.FIPS, row.SCC, row.POLID, row.date)] = row.hourly_values
    end
    return cache
end

"""
    temporal_allocate(emissions::DataFrame, profiles::DataFrame,
        xref::DataFrame, episode_start::DateTime, episode_end::DateTime;
        timezone_offset::Int=0,
        timezone_map::Dict{String,Int}=Dict{String,Int}(),
        day_specific::DataFrame=DataFrame(),
        hour_specific::DataFrame=DataFrame()) -> DataFrame

Convert annual emissions to hourly emissions using temporal profiles.

Applies monthly, day-of-week, and diurnal temporal factors to distribute
annual emissions to hourly values for the specified episode period.

Supports SMOKE-style priority hierarchy for temporal overrides:
1. Hour-specific data → use directly (bypass all profiles)
2. Day-specific data → use with diurnal profile (bypass monthly/weekly)
3. Annual with profiles → standard monthly × weekly × diurnal allocation

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
- `day_specific::DataFrame=DataFrame()`: Day-specific emissions from [`read_day_specific`](@ref).
  When matched, overrides monthly/weekly profile allocation.
- `hour_specific::DataFrame=DataFrame()`: Hour-specific emissions from [`read_hour_specific`](@ref).
  When matched, overrides all profile-based allocation.

# Returns
A `DataFrame` with columns: `:FIPS`, `:SCC`, `:POLID`, `:hour` (DateTime),
`:emission_rate` (Float64), plus any location columns from the input.

The temporal allocation follows SMOKE conventions:
- Monthly factor: fraction of annual total for the given month (sums to 1.0)
- Weekly factor: relative weight for the day of week (sums to 7.0)
- Diurnal factor: fraction of daily total for the given hour (sums to 1.0)
- Hourly rate = ANN_VALUE × (monthly_factor × 12) × weekly_factor × (diurnal_factor × 24)
  This converts from fraction-based profiles to rate multipliers:
  - monthly_factor × 12: converts from "fraction of annual" to rate modifier (uniform → 1.0)
  - weekly_factor: day-of-week relative weight (uniform → 1.0, double-weight → 2.0)
  - diurnal_factor × 24: converts from "fraction of daily" to rate modifier (uniform → 1.0)
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
        timezone_map::Dict{String, Int} = Dict{String, Int}(),
        day_specific::DataFrame = DataFrame(),
        hour_specific::DataFrame = DataFrame()
    )
    # Pre-compute hours in the episode
    n_hours = Int(Dates.value(episode_end - episode_start)) ÷ 3_600_000
    hours = [episode_start + Hour(h) for h in 0:(n_hours - 1)]

    # Detect extra columns to carry through
    extra_cols = Symbol[]
    for col in [:LONGITUDE, :LATITUDE, :COUNTRY, :Surrogate]
        if hasproperty(emissions, col)
            push!(extra_cols, col)
        end
    end

    if isempty(hours) || nrow(emissions) == 0
        result = DataFrame(
            FIPS = String[],
            SCC = String[],
            POLID = String[],
            hour = DateTime[],
            emission_rate = Float64[]
        )
        for col in extra_cols
            result[!, col] = Any[]
        end
        return result
    end

    n_sources = nrow(emissions)
    n_total = n_sources * length(hours)

    # Build caches for O(1) lookups
    profile_cache = _build_profile_cache(profiles)
    xref_cache = _build_xref_cache(xref)

    # Build day/hour-specific caches if provided
    has_day_specific = nrow(day_specific) > 0
    has_hour_specific = nrow(hour_specific) > 0
    day_cache = has_day_specific ? _build_day_specific_cache(day_specific) : _DaySpecificCache()
    hour_cache = has_hour_specific ? _build_hour_specific_cache(hour_specific) : _HourSpecificCache()

    # Pre-allocate output column vectors
    out_fips = Vector{String}(undef, n_total)
    out_scc = Vector{String}(undef, n_total)
    out_polid = Vector{String}(undef, n_total)
    out_hour = Vector{DateTime}(undef, n_total)
    out_rate = Vector{Float64}(undef, n_total)
    out_extra = Dict{Symbol, Vector{Any}}(col => Vector{Any}(undef, n_total) for col in extra_cols)

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

        # Cache extra column values for this source
        extra_vals = Dict{Symbol, Any}(col => erow[col] for col in extra_cols)

        for hr in hours
            local_hr = hr + Hour(tz_offset)
            month_idx = Dates.month(local_hr)
            dow = Dates.dayofweek(local_hr)  # 1=Monday, 7=Sunday
            hour_idx = Dates.hour(local_hr) + 1  # 1-indexed
            local_date = Date(local_hr)

            hourly_rate = 0.0

            # Priority 1: Hour-specific override
            if has_hour_specific
                hour_key = (fips, scc, polid, local_date)
                hour_vals = get(hour_cache, hour_key, nothing)
                if hour_vals !== nothing && hour_idx <= length(hour_vals)
                    hourly_rate = hour_vals[hour_idx]
                    idx += 1
                    out_fips[idx] = fips
                    out_scc[idx] = scc
                    out_polid[idx] = polid
                    out_hour[idx] = hr
                    out_rate[idx] = hourly_rate
                    for col in extra_cols
                        out_extra[col][idx] = extra_vals[col]
                    end
                    continue
                end
            end

            # Priority 2: Day-specific override (still uses diurnal profile)
            if has_day_specific
                day_key = (fips, scc, polid, local_date)
                day_val = get(day_cache, day_key, nothing)
                if day_val !== nothing
                    diurnal_factors = _lookup_diurnal_profile(profile_cache, profile_ids.diurnal_id, dow)
                    df = hour_idx <= length(diurnal_factors) ? diurnal_factors[hour_idx] : 1.0 / 24.0
                    # Day value distributed across hours using diurnal profile
                    hourly_rate = day_val * df * 24.0
                    idx += 1
                    out_fips[idx] = fips
                    out_scc[idx] = scc
                    out_polid[idx] = polid
                    out_hour[idx] = hr
                    out_rate[idx] = hourly_rate
                    for col in extra_cols
                        out_extra[col][idx] = extra_vals[col]
                    end
                    continue
                end
            end

            # Priority 3: Standard annual allocation with profiles
            diurnal_factors = _lookup_diurnal_profile(profile_cache, profile_ids.diurnal_id, dow)

            mf = month_idx <= length(monthly_factors) ? monthly_factors[month_idx] : 1.0 / 12.0
            wf = dow <= length(weekly_factors) ? weekly_factors[dow] : 1.0
            df = hour_idx <= length(diurnal_factors) ? diurnal_factors[hour_idx] : 1.0 / 24.0

            # hourly_rate = annual_rate * (monthly_frac*12) * weekly_weight * (diurnal_frac * 24)
            hourly_rate = ann_value * (mf * 12.0) * wf * (df * 24.0)

            idx += 1
            out_fips[idx] = fips
            out_scc[idx] = scc
            out_polid[idx] = polid
            out_hour[idx] = hr
            out_rate[idx] = hourly_rate
            for col in extra_cols
                out_extra[col][idx] = extra_vals[col]
            end
        end
    end

    result = DataFrame(
        FIPS = out_fips,
        SCC = out_scc,
        POLID = out_polid,
        hour = out_hour,
        emission_rate = out_rate,
    )
    for col in extra_cols
        result[!, col] = out_extra[col]
    end
    return result
end
