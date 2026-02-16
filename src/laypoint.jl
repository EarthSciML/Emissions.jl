export LayerConfig, MetProfile, compute_layer_fractions,
    allocate_point_to_layers, laypoint

"""
    LayerConfig

Configuration for vertical atmospheric layers.

# Fields
- `layer_heights::Vector{Float64}`: Staggered layer interface heights in meters
  (n+1 values for n layers). First value is surface (usually 0.0), last is model top.
- `n_layers::Int`: Number of atmospheric layers.
"""
struct LayerConfig
    layer_heights::Vector{Float64}
    n_layers::Int

    function LayerConfig(heights::Vector{Float64})
        n = length(heights) - 1
        n >= 1 || throw(ArgumentError("Need at least 2 height values for 1 layer"))
        return new(heights, n)
    end
end

"""
    MetProfile

Meteorological profile data for vertical layer allocation.

All vectors should have length equal to the number of layers.

# Fields
- `temperature::Vector{Float64}`: Temperature at each layer center (K)
- `wind_speed::Vector{Float64}`: Wind speed at each layer center (m/s)
- `stability_class::Vector{Float64}`: Stability class (0=unstable, >0.5=stable)
- `stability_param::Vector{Float64}`: Stability parameter s1 for plume rise
"""
struct MetProfile
    temperature::Vector{Float64}
    wind_speed::Vector{Float64}
    stability_class::Vector{Float64}
    stability_param::Vector{Float64}
end

"""
    compute_layer_fractions(plume_bottom::Float64, plume_top::Float64,
        config::LayerConfig) -> Vector{Float64}

Compute the fraction of a plume within each atmospheric layer using
pressure-weighted allocation.

The plume is assumed to be uniformly distributed between `plume_bottom`
and `plume_top`. The fraction in each layer is proportional to the overlap
height between the plume and that layer.

Fractions sum to 1.0. If the plume extends above the model top, the
excess is allocated to the top layer.

# Arguments
- `plume_bottom`: Bottom of the plume (m), typically the stack height.
- `plume_top`: Top of the plume (m), stack height + plume rise.
- `config`: Layer configuration.

# Returns
A `Vector{Float64}` of length `config.n_layers` with layer fractions summing to 1.0.
"""
function compute_layer_fractions(
        plume_bottom::Float64, plume_top::Float64,
        config::LayerConfig,
    )
    fracs = zeros(Float64, config.n_layers)

    # Clamp to valid range
    plume_bottom = max(plume_bottom, config.layer_heights[1])
    plume_top = max(plume_top, plume_bottom)

    plume_depth = plume_top - plume_bottom

    if plume_depth <= 0.0
        # Point source at a single height - find the containing layer
        for k in 1:config.n_layers
            if plume_bottom >= config.layer_heights[k] && plume_bottom < config.layer_heights[k + 1]
                fracs[k] = 1.0
                return fracs
            end
        end
        # Above model top -> put in top layer
        fracs[config.n_layers] = 1.0
        return fracs
    end

    # Distribute plume across layers proportionally
    for k in 1:config.n_layers
        layer_bot = config.layer_heights[k]
        layer_top = config.layer_heights[k + 1]

        overlap_bot = max(plume_bottom, layer_bot)
        overlap_top = min(plume_top, layer_top)

        if overlap_top > overlap_bot
            fracs[k] = (overlap_top - overlap_bot) / plume_depth
        end
    end

    # If plume extends above model top, assign remaining to top layer
    total = sum(fracs)
    if total < 1.0 && plume_top > config.layer_heights[end]
        fracs[config.n_layers] += (1.0 - total)
    end

    # Normalize to ensure exact sum of 1.0
    s = sum(fracs)
    if s > 0.0
        fracs ./= s
    end

    return fracs
end

"""
    allocate_point_to_layers(height, diam, temp, vel, met::MetProfile,
        config::LayerConfig; spread::Float64=0.5) -> Vector{Float64}

Compute layer fractions for a single point source using ASME plume rise.

Calculates plume rise using the existing [`ASME`](@ref) algorithm with
the provided meteorological profile, then distributes the plume across
layers using [`compute_layer_fractions`](@ref).

The plume is distributed over a range from `stack_height` to
`stack_height + plume_rise + spread * plume_rise`.

# Arguments
- `height`: Stack height (m)
- `diam`: Stack diameter (m)
- `temp`: Stack exit temperature (K)
- `vel`: Stack exit velocity (m/s)
- `met::MetProfile`: Meteorological profile for the source location.
- `config::LayerConfig`: Layer configuration.
- `spread::Float64=0.5`: Plume spread factor (fraction of plume rise added above center).

# Returns
Layer fractions as `Vector{Float64}` summing to 1.0.
"""
function allocate_point_to_layers(
        height::Real, diam::Real, temp::Real, vel::Real,
        met::MetProfile, config::LayerConfig;
        spread::Float64 = 0.5,
    )
    height = Float64(height)
    diam = Float64(diam)
    temp = Float64(temp)
    vel = Float64(vel)

    # Use ASME plume rise calculation
    plume_layer = 0
    plume_height = height
    try
        plume_layer, plume_height = ASME(
            height, diam, temp, vel,
            config.layer_heights,
            met.temperature, met.wind_speed,
            met.stability_class, met.stability_param,
        )
    catch e
        if e == ErrAboveModelTop
            # Plume above model top - put everything in top layer
            fracs = zeros(Float64, config.n_layers)
            fracs[config.n_layers] = 1.0
            return fracs
        else
            rethrow(e)
        end
    end

    delta_h = plume_height - height
    plume_bottom = height
    plume_top = plume_height + spread * max(delta_h, 0.0)

    return compute_layer_fractions(plume_bottom, plume_top, config)
end

"""
    laypoint(point_emissions::DataFrame, met_profiles::Dict{String, MetProfile},
        config::LayerConfig) -> DataFrame

Allocate point source emissions to vertical atmospheric layers.

For each point source, computes plume rise using the ASME algorithm and
distributes emissions across layers. Sources without meteorological data
are placed in the surface layer.

# Arguments
- `point_emissions::DataFrame`: Must have columns `:STKHGT`, `:STKDIAM`,
  `:STKTEMP`, `:STKVEL`, `:FIPS`. May have `:LONGITUDE`, `:LATITUDE`.
- `met_profiles::Dict{String, MetProfile}`: Meteorological profiles keyed by
  location key (see [`location_key`](@ref)). A key of `"default"` is used as fallback.
- `config::LayerConfig`: Layer configuration.

# Returns
A `DataFrame` with all original columns plus:
- `:layer` (Int): Layer index (1-based)
- `:layer_fraction` (Float64): Fraction of emissions in this layer

Each input row may produce multiple output rows (one per layer with non-zero fraction).
"""
function laypoint(
        point_emissions::DataFrame,
        met_profiles::Dict{String, MetProfile},
        config::LayerConfig,
    )
    has_lon = hasproperty(point_emissions, :LONGITUDE)
    has_lat = hasproperty(point_emissions, :LATITUDE)

    # Collect output columns
    all_cols = names(point_emissions)
    out_data = Dict{String, Vector{Any}}(col => Any[] for col in all_cols)
    out_layer = Int[]
    out_frac = Float64[]

    for erow in eachrow(point_emissions)
        height = Float64(erow.STKHGT)
        diam = Float64(erow.STKDIAM)
        temp = Float64(erow.STKTEMP)
        vel = Float64(erow.STKVEL)

        # Build location key for met profile lookup
        fips = string(erow.FIPS)
        lon = has_lon ? erow.LONGITUDE : missing
        lat = has_lat ? erow.LATITUDE : missing
        key = location_key(fips, lon, lat)

        # Get met profile (try source-specific, then default)
        met = get(met_profiles, key, get(met_profiles, "default", nothing))

        if met !== nothing
            fracs = allocate_point_to_layers(height, diam, temp, vel, met, config)
        else
            # No met data - put in surface layer
            fracs = zeros(Float64, config.n_layers)
            # Find layer containing stack
            for k in 1:config.n_layers
                if height >= config.layer_heights[k] && height < config.layer_heights[k + 1]
                    fracs[k] = 1.0
                    break
                end
            end
            if sum(fracs) == 0.0
                fracs[1] = 1.0
            end
        end

        # Expand into output rows (one per non-zero layer fraction)
        for k in 1:config.n_layers
            if fracs[k] > 0.0
                for col in all_cols
                    push!(out_data[col], erow[Symbol(col)])
                end
                push!(out_layer, k)
                push!(out_frac, fracs[k])
            end
        end
    end

    result = DataFrame()
    for col in all_cols
        result[!, Symbol(col)] = out_data[col]
    end
    result[!, :layer] = out_layer
    result[!, :layer_fraction] = out_frac

    return result
end
