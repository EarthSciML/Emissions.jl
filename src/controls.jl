export ControlSpec, read_growth_factors, read_control_factors, apply_controls

"""
    ControlSpec

Specification for an emissions control or growth factor.

# Fields
- `region::String`: Country or region code
- `fips::String`: 5-digit FIPS code ("00000" for national default)
- `scc::String`: Source Classification Code ("0000000000" for any SCC)
- `pollutant::String`: Pollutant name ("" for all pollutants)
- `control_type::Symbol`: `:growth`, `:multiplicative`, or `:reactivity`
- `factor::Float64`: Multiplicative factor to apply to emissions
- `base_year::Int`: Base year of the inventory
- `target_year::Int`: Target year for projection
"""
struct ControlSpec
    region::String
    fips::String
    scc::String
    pollutant::String
    control_type::Symbol
    factor::Float64
    base_year::Int
    target_year::Int
end

"""
    read_growth_factors(filepath::AbstractString) -> Vector{ControlSpec}

Read a SMOKE-format growth factor file.

Lines starting with `#` are skipped. Each data line has semicolon-delimited fields:
`FIPS;SCC;pollutant;base_year;target_year;growth_factor`

Growth factors are applied multiplicatively: `new_value = old_value * growth_factor`.

# Returns
A `Vector{ControlSpec}` with `control_type = :growth`.
"""
function read_growth_factors(filepath::AbstractString)
    lines = readlines(filepath)
    controls = ControlSpec[]
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        # Remove inline comments after '!'
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ';')
        if length(parts) >= 6
            fips_raw = strip(parts[1])
            scc = strip(parts[2])
            pollutant = uppercase(strip(parts[3]))
            base_year = parse(Int, strip(parts[4]))
            target_year = parse(Int, strip(parts[5]))
            factor = parse(Float64, strip(parts[6]))

            # Normalize FIPS
            if length(fips_raw) == 6
                fips = lpad(fips_raw[2:end], 5, '0')
            else
                fips = lpad(fips_raw, 5, '0')
            end

            # Normalize SCC
            scc = lpad(scc, 10, '0')

            push!(
                controls, ControlSpec(
                    "", fips, scc, pollutant, :growth, factor, base_year, target_year
                )
            )
        end
    end
    return controls
end

"""
    read_control_factors(filepath::AbstractString) -> Vector{ControlSpec}

Read a SMOKE-format multiplicative control factor file.

Lines starting with `#` are skipped. Each data line has semicolon-delimited fields:
`FIPS;SCC;pollutant;efficiency;effectiveness;penetration`

The multiplicative control factor is computed as:
`factor = 1 - (efficiency/100 * effectiveness/100 * penetration/100)`

# Returns
A `Vector{ControlSpec}` with `control_type = :multiplicative`.
"""
function read_control_factors(filepath::AbstractString)
    lines = readlines(filepath)
    controls = ControlSpec[]
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        # Remove inline comments after '!'
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ';')
        if length(parts) >= 6
            fips_raw = strip(parts[1])
            scc = strip(parts[2])
            pollutant = uppercase(strip(parts[3]))
            efficiency = parse(Float64, strip(parts[4]))
            effectiveness = parse(Float64, strip(parts[5]))
            penetration = parse(Float64, strip(parts[6]))

            # Normalize FIPS
            if length(fips_raw) == 6
                fips = lpad(fips_raw[2:end], 5, '0')
            else
                fips = lpad(fips_raw, 5, '0')
            end

            # Normalize SCC
            scc = lpad(scc, 10, '0')

            factor = 1.0 - (efficiency / 100.0 * effectiveness / 100.0 * penetration / 100.0)

            push!(
                controls, ControlSpec(
                    "", fips, scc, pollutant, :multiplicative, factor, 0, 0
                )
            )
        end
    end
    return controls
end

"""
    _match_control(controls::Vector{ControlSpec}, fips::AbstractString,
        scc::AbstractString, pollutant::AbstractString,
        control_type::Symbol) -> Union{ControlSpec, Nothing}

Match an emissions source to a control specification using 7-level hierarchical matching:
1. Exact FIPS + SCC + pollutant
2. Exact FIPS + SCC (any pollutant)
3. Exact FIPS + pollutant (any SCC)
4. Exact FIPS (any SCC, any pollutant)
5. National default (FIPS="00000") + SCC + pollutant
6. National default + pollutant (any SCC)
7. National default (any SCC, any pollutant)

Only controls matching the specified `control_type` are considered.
"""
function _match_control(
        controls::Vector{ControlSpec},
        fips::AbstractString,
        scc::AbstractString,
        pollutant::AbstractString,
        control_type::Symbol,
    )
    pollutant_up = uppercase(pollutant)

    # Filter to matching type
    typed = [c for c in controls if c.control_type == control_type]

    # Level 1: FIPS + SCC + pollutant
    for c in typed
        if c.fips == fips && c.scc == scc && c.pollutant == pollutant_up
            return c
        end
    end

    # Level 2: FIPS + SCC (any pollutant)
    for c in typed
        if c.fips == fips && c.scc == scc && c.pollutant == ""
            return c
        end
    end

    # Level 3: FIPS + pollutant (any SCC)
    for c in typed
        if c.fips == fips && c.scc == "0000000000" && c.pollutant == pollutant_up
            return c
        end
    end

    # Level 4: FIPS only (any SCC, any pollutant)
    for c in typed
        if c.fips == fips && c.scc == "0000000000" && c.pollutant == ""
            return c
        end
    end

    # Level 5: National default + SCC + pollutant
    for c in typed
        if c.fips == "00000" && c.scc == scc && c.pollutant == pollutant_up
            return c
        end
    end

    # Level 6: National default + pollutant (any SCC)
    for c in typed
        if c.fips == "00000" && c.scc == "0000000000" && c.pollutant == pollutant_up
            return c
        end
    end

    # Level 7: National default (any SCC, any pollutant)
    for c in typed
        if c.fips == "00000" && c.scc == "0000000000" && c.pollutant == ""
            return c
        end
    end

    return nothing
end

"""
    apply_controls(emissions::DataFrame, controls::Vector{ControlSpec}) -> DataFrame

Apply emissions controls to scale `ANN_VALUE`.

Growth factors are applied first (multiplicative), then multiplicative controls.
The combined factor is: `new_ANN_VALUE = ANN_VALUE * growth_factor * control_factor`.

# Arguments
- `emissions::DataFrame`: Must have columns `:FIPS`, `:SCC`, `:POLID`, `:ANN_VALUE`.
- `controls::Vector{ControlSpec}`: Control specifications from [`read_growth_factors`](@ref)
  and/or [`read_control_factors`](@ref).

# Returns
A copy of `emissions` with `ANN_VALUE` adjusted by matching control factors.
"""
function apply_controls(emissions::DataFrame, controls::Vector{ControlSpec})
    result = copy(emissions)

    for i in 1:nrow(result)
        fips = string(result[i, :FIPS])
        scc = lpad(string(result[i, :SCC]), 10, '0')
        polid = uppercase(string(result[i, :POLID]))

        combined_factor = 1.0

        # Apply growth factors first
        growth = _match_control(controls, fips, scc, polid, :growth)
        if growth !== nothing
            combined_factor *= growth.factor
        end

        # Apply multiplicative controls
        mult = _match_control(controls, fips, scc, polid, :multiplicative)
        if mult !== nothing
            combined_factor *= mult.factor
        end

        result[i, :ANN_VALUE] = result[i, :ANN_VALUE] * combined_factor
    end

    return result
end
