export ElevationCriteria, DEFAULT_ELEVATION_CRITERIA,
    analytical_plume_rise, classify_point_sources, group_stacks

"""
    ElevationCriteria

Thresholds for classifying point sources as surface, elevated, or
Plume-in-Grid (PinG).

All values in SI units (meters, m/s, Kelvin).

# Fields
- `min_stack_height::Float64`: Minimum stack height for elevated (m)
- `min_exit_velocity::Float64`: Minimum exit velocity for elevated (m/s)
- `min_exit_temperature::Float64`: Minimum exit temperature for elevated (K)
- `min_flow_rate::Float64`: Minimum volumetric flow rate for elevated (m³/s)
- `min_plume_rise::Float64`: Minimum analytical plume rise for elevated (m)
- `ping_stack_height::Float64`: Minimum stack height for PinG (m)
- `ping_emissions_threshold::Float64`: Minimum emissions rate for PinG (kg/s)
"""
struct ElevationCriteria
    min_stack_height::Float64
    min_exit_velocity::Float64
    min_exit_temperature::Float64
    min_flow_rate::Float64
    min_plume_rise::Float64
    ping_stack_height::Float64
    ping_emissions_threshold::Float64
end

"""
    DEFAULT_ELEVATION_CRITERIA

SMOKE default thresholds for elevated source classification.
Converted from imperial to SI:
- 50 ft = 15.24 m (stack height)
- 10 ft/s = 3.048 m/s (exit velocity)
- 200°F = 366.48 K (exit temperature)
- 100 ft³/s = 2.8317 m³/s (flow rate)
- 100 ft = 30.48 m (plume rise and PinG height)
"""
const DEFAULT_ELEVATION_CRITERIA = ElevationCriteria(
    15.24,   # min_stack_height (50 ft)
    3.048,   # min_exit_velocity (10 ft/s)
    366.48,  # min_exit_temperature (200°F)
    2.8317,  # min_flow_rate (100 ft³/s)
    30.48,   # min_plume_rise (100 ft)
    30.48,   # ping_stack_height (100 ft)
    0.0,     # ping_emissions_threshold (no minimum)
)

"""
    analytical_plume_rise(height, diam, temp, vel;
        ambient_temp=293.15, ambient_wind=3.0) -> Float64

Compute Briggs analytical plume rise for screening purposes.

Uses default meteorological values (ambient temperature 293.15 K,
wind speed 3.0 m/s) as in SMOKE's Elevpoint program.

# Arguments
- `height`: Stack height (m)
- `diam`: Stack diameter (m)
- `temp`: Stack exit temperature (K)
- `vel`: Stack exit velocity (m/s)
- `ambient_temp`: Ambient temperature (K), default 293.15
- `ambient_wind`: Ambient wind speed (m/s), default 3.0

# Returns
Analytical plume rise in meters.
"""
function analytical_plume_rise(
        height::Real, diam::Real, temp::Real, vel::Real;
        ambient_temp::Real = 293.15, ambient_wind::Real = 3.0,
    )
    height = Float64(height)
    diam = Float64(diam)
    temp = Float64(temp)
    vel = Float64(vel)
    ambient_temp = Float64(ambient_temp)
    ambient_wind = Float64(ambient_wind)

    if ambient_wind <= 0.0
        ambient_wind = 1.0
    end

    temp_diff = temp - ambient_temp

    if temp_diff < 50.0 && vel > ambient_wind && vel > 10.0
        # Momentum-dominated plume rise
        deltaH = diam * (vel^1.4) / (ambient_wind^1.4)
    else
        # Buoyancy-dominated plume rise
        F = g * max(temp_diff, 0.0) / max(temp, 1.0) * vel * (diam / 2.0)^2

        if F > 0.0
            # Unstable conditions (default for screening)
            deltaH = 7.4 * (F * height^2.0)^(1.0 / 3.0) / ambient_wind
        else
            deltaH = 0.0
        end
    end

    return max(deltaH, 0.0)
end

"""
    classify_point_sources(point_emissions::DataFrame;
        criteria::ElevationCriteria=DEFAULT_ELEVATION_CRITERIA) -> DataFrame

Classify point sources as "surface", "elevated", or "ping" (Plume-in-Grid).

A source is classified as "elevated" if it meets ANY of these criteria:
- Stack height ≥ `min_stack_height`
- Exit velocity ≥ `min_exit_velocity`
- Exit temperature ≥ `min_exit_temperature`
- Flow rate ≥ `min_flow_rate`
- Analytical plume rise ≥ `min_plume_rise`

A source is further classified as "ping" if it meets ALL PinG criteria:
- Stack height ≥ `ping_stack_height`
- Emissions rate ≥ `ping_emissions_threshold` (if threshold > 0)

# Arguments
- `point_emissions::DataFrame`: Must have columns `:STKHGT` (m), `:STKDIAM` (m),
  `:STKTEMP` (K), `:STKVEL` (m/s). May have `:ANN_VALUE` for PinG classification.
- `criteria`: Elevation criteria thresholds.

# Returns
A copy of `point_emissions` with added columns:
- `:source_class` (String): "surface", "elevated", or "ping"
- `:analytical_plume_rise` (Float64): Computed analytical plume rise (m)
"""
function classify_point_sources(
        point_emissions::DataFrame;
        criteria::ElevationCriteria = DEFAULT_ELEVATION_CRITERIA,
    )
    result = copy(point_emissions)
    n = nrow(result)

    source_class = fill("surface", n)
    plume_rises = zeros(Float64, n)

    has_ann = hasproperty(result, :ANN_VALUE)

    for i in 1:n
        height = Float64(result[i, :STKHGT])
        diam = Float64(result[i, :STKDIAM])
        temp = Float64(result[i, :STKTEMP])
        vel = Float64(result[i, :STKVEL])

        # Compute flow rate (m³/s)
        flow_rate = vel * π * (diam / 2.0)^2

        # Compute analytical plume rise
        pr = analytical_plume_rise(height, diam, temp, vel)
        plume_rises[i] = pr

        # Check elevated criteria (ANY threshold)
        is_elevated = (
            height >= criteria.min_stack_height ||
                vel >= criteria.min_exit_velocity ||
                temp >= criteria.min_exit_temperature ||
                flow_rate >= criteria.min_flow_rate ||
                pr >= criteria.min_plume_rise
        )

        if is_elevated
            source_class[i] = "elevated"

            # Check PinG criteria (ALL thresholds)
            is_ping = height >= criteria.ping_stack_height
            if is_ping && criteria.ping_emissions_threshold > 0.0 && has_ann
                ann_val = Float64(result[i, :ANN_VALUE] isa Number ? result[i, :ANN_VALUE] : ustrip(result[i, :ANN_VALUE]))
                is_ping = ann_val >= criteria.ping_emissions_threshold
            end

            if is_ping
                source_class[i] = "ping"
            end
        end
    end

    result[!, :source_class] = source_class
    result[!, :analytical_plume_rise] = plume_rises
    return result
end

"""
    group_stacks(point_emissions::DataFrame;
        height_bin::Float64=10.0, temp_bin::Float64=50.0) -> DataFrame

Group similar stacks at the same facility based on stack parameters.

Stacks at the same facility (same FIPS + facility ID or same coordinates)
are grouped if their stack height and temperature are within the specified
bin widths.

# Arguments
- `point_emissions::DataFrame`: Must have `:STKHGT`, `:STKTEMP`, `:FIPS`.
  May have `:LONGITUDE`, `:LATITUDE` for location-based grouping.
- `height_bin`: Stack height bin width in meters (default 10.0).
- `temp_bin`: Temperature bin width in Kelvin (default 50.0).

# Returns
A copy of `point_emissions` with added column `:stack_group` (Int).
"""
function group_stacks(
        point_emissions::DataFrame;
        height_bin::Float64 = 10.0, temp_bin::Float64 = 50.0,
    )
    result = copy(point_emissions)
    n = nrow(result)
    stack_group = zeros(Int, n)
    group_counter = 0

    has_lon = hasproperty(result, :LONGITUDE)
    has_lat = hasproperty(result, :LATITUDE)

    for i in 1:n
        if stack_group[i] != 0
            continue
        end
        group_counter += 1
        stack_group[i] = group_counter

        fips_i = string(result[i, :FIPS])
        height_i = Float64(result[i, :STKHGT])
        temp_i = Float64(result[i, :STKTEMP])

        lon_i = has_lon ? result[i, :LONGITUDE] : missing
        lat_i = has_lat ? result[i, :LATITUDE] : missing

        for j in (i + 1):n
            if stack_group[j] != 0
                continue
            end

            fips_j = string(result[j, :FIPS])
            lon_j = has_lon ? result[j, :LONGITUDE] : missing
            lat_j = has_lat ? result[j, :LATITUDE] : missing

            # Same location check
            same_location = false
            if !ismissing(lon_i) && !ismissing(lat_i) && !ismissing(lon_j) && !ismissing(lat_j)
                same_location = (lon_i == lon_j && lat_i == lat_j && fips_i == fips_j)
            else
                same_location = (fips_i == fips_j)
            end

            if !same_location
                continue
            end

            height_j = Float64(result[j, :STKHGT])
            temp_j = Float64(result[j, :STKTEMP])

            if abs(height_i - height_j) <= height_bin && abs(temp_i - temp_j) <= temp_bin
                stack_group[j] = group_counter
            end
        end
    end

    result[!, :stack_group] = stack_group
    return result
end
