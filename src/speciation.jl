export read_gspro, read_gsref, build_speciation_matrix, speciate_emissions

"""
    read_gspro(filepath::AbstractString) -> DataFrame

Read a SMOKE-format GSPRO speciation profile file.

Lines starting with `#` are skipped. Each data line has semicolon-delimited fields:
`profile_code;pollutant_id;species_id;split_factor;divisor;mass_fraction`

# Returns
A `DataFrame` with columns:
- `:profile_code` (String): Speciation profile identifier
- `:pollutant_id` (String): Inventory pollutant name (e.g., "VOC", "NOX")
- `:species_id` (String): Model species name (e.g., "NO", "NO2", "FORM")
- `:split_factor` (Float64): Mole-based split factor (numerator)
- `:divisor` (Float64): Mole-based divisor (denominator)
- `:mass_fraction` (Float64): Mass-based speciation fraction
"""
function read_gspro(filepath::AbstractString)
    lines = readlines(filepath)
    records = DataFrame(
        profile_code = String[],
        pollutant_id = String[],
        species_id = String[],
        split_factor = Float64[],
        divisor = Float64[],
        mass_fraction = Float64[],
    )
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
            profile_code = strip(parts[1])
            pollutant_id = uppercase(strip(parts[2]))
            species_id = strip(parts[3])
            split_factor = parse(Float64, strip(parts[4]))
            divisor = parse(Float64, strip(parts[5]))
            mass_fraction = parse(Float64, strip(parts[6]))
            push!(
                records, (
                    profile_code = profile_code,
                    pollutant_id = pollutant_id,
                    species_id = species_id,
                    split_factor = split_factor,
                    divisor = divisor,
                    mass_fraction = mass_fraction,
                )
            )
        end
    end
    return records
end

"""
    read_gsref(filepath::AbstractString) -> DataFrame

Read a SMOKE-format GSREF speciation cross-reference file.

Maps emission sources (by FIPS and SCC) to speciation profiles for each pollutant.
Lines starting with `#` are skipped. Each data line has semicolon-delimited fields:
`FIPS;SCC;pollutant_id;profile_code`

FIPS codes are normalized: 6-digit codes have the country digit stripped and are
left-padded to 5 digits. `"00000"` serves as the national default.

# Returns
A `DataFrame` with columns:
- `:FIPS` (String): 5-digit FIPS code ("00000" for national default)
- `:SCC` (String): Source Classification Code
- `:pollutant_id` (String): Inventory pollutant name
- `:profile_code` (String): Speciation profile code (matches GSPRO)
"""
function read_gsref(filepath::AbstractString)
    lines = readlines(filepath)
    records = DataFrame(
        FIPS = String[],
        SCC = String[],
        pollutant_id = String[],
        profile_code = String[],
    )
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
        if length(parts) >= 4
            fips_raw = strip(parts[1])
            scc = strip(parts[2])
            pollutant_id = uppercase(strip(parts[3]))
            profile_code = strip(parts[4])

            # Normalize FIPS
            if length(fips_raw) == 6
                fips = lpad(fips_raw[2:end], 5, '0')
            else
                fips = lpad(fips_raw, 5, '0')
            end

            # Normalize SCC to 10 digits
            scc = lpad(scc, 10, '0')

            push!(
                records, (
                    FIPS = fips,
                    SCC = scc,
                    pollutant_id = pollutant_id,
                    profile_code = profile_code,
                )
            )
        end
    end
    return records
end

"""
    _match_speciation_profile(gsref::DataFrame, fips::AbstractString,
        scc::AbstractString, pollutant::AbstractString) -> Union{String, Nothing}

Match an emissions source to a speciation profile using hierarchical matching:
1. Exact FIPS + SCC + pollutant match
2. National default (FIPS="00000") + SCC + pollutant match
3. Pollutant-only match (any FIPS, any SCC where SCC is "0000000000")
4. `nothing` if no match found
"""
function _match_speciation_profile(
        gsref::DataFrame,
        fips::AbstractString,
        scc::AbstractString,
        pollutant::AbstractString,
    )
    pollutant_up = uppercase(pollutant)

    # Level 1: Exact FIPS + SCC + pollutant
    for row in eachrow(gsref)
        if row.FIPS == fips && row.SCC == scc && row.pollutant_id == pollutant_up
            return row.profile_code
        end
    end

    # Level 2: State-level FIPS (first 2 digits + "000") + SCC + pollutant
    if length(fips) >= 2
        state_fips = fips[1:2] * "000"
        for row in eachrow(gsref)
            if row.FIPS == state_fips && row.SCC == scc && row.pollutant_id == pollutant_up
                return row.profile_code
            end
        end
    end

    # Level 3: National default FIPS + SCC + pollutant
    for row in eachrow(gsref)
        if row.FIPS == "00000" && row.SCC == scc && row.pollutant_id == pollutant_up
            return row.profile_code
        end
    end

    # Level 4: Pollutant-only (any FIPS, SCC="0000000000")
    for row in eachrow(gsref)
        if row.SCC == "0000000000" && row.pollutant_id == pollutant_up
            return row.profile_code
        end
    end

    return nothing
end

"""
    build_speciation_matrix(emissions::DataFrame, gspro::DataFrame, gsref::DataFrame;
        basis::Symbol=:mass) -> (Matrix{Float64}, Vector{String}, Vector{Int})

Build a speciation matrix mapping inventory sources to model species.

# Arguments
- `emissions::DataFrame`: Must have columns `:FIPS`, `:SCC`, `:POLID`.
- `gspro::DataFrame`: Speciation profiles from [`read_gspro`](@ref).
- `gsref::DataFrame`: Cross-reference from [`read_gsref`](@ref).
- `basis::Symbol=:mass`: `:mass` uses `mass_fraction`, `:mole` uses `split_factor/divisor`.

# Returns
A tuple of:
- `matrix::Matrix{Float64}`: (n_sources × n_species) speciation factors
- `species_names::Vector{String}`: Species names for each column
- `source_indices::Vector{Int}`: Row index in `emissions` for each matrix row
"""
function build_speciation_matrix(
        emissions::DataFrame,
        gspro::DataFrame,
        gsref::DataFrame;
        basis::Symbol = :mass,
    )
    # Collect all unique species from GSPRO
    all_species = unique(gspro.species_id)
    sort!(all_species)
    species_idx = Dict(sp => i for (i, sp) in enumerate(all_species))
    n_species = length(all_species)

    source_indices = Int[]
    rows_data = Vector{Float64}[]

    for (i, erow) in enumerate(eachrow(emissions))
        fips = string(erow.FIPS)
        scc = lpad(string(erow.SCC), 10, '0')
        polid = uppercase(string(erow.POLID))

        profile_code = _match_speciation_profile(gsref, fips, scc, polid)
        if profile_code === nothing
            continue
        end

        # Look up profile entries in GSPRO
        row_vec = zeros(Float64, n_species)
        found = false
        for prow in eachrow(gspro)
            if prow.profile_code == profile_code && prow.pollutant_id == polid
                sp_idx = get(species_idx, prow.species_id, 0)
                if sp_idx > 0
                    if basis == :mole
                        factor = prow.divisor != 0.0 ? prow.split_factor / prow.divisor : 0.0
                    else
                        factor = prow.mass_fraction
                    end
                    row_vec[sp_idx] = factor
                    found = true
                end
            end
        end

        if found
            push!(source_indices, i)
            push!(rows_data, row_vec)
        end
    end

    if isempty(rows_data)
        return zeros(Float64, 0, n_species), all_species, Int[]
    end

    matrix = reduce(vcat, [r' for r in rows_data])
    return matrix, all_species, source_indices
end

"""
    speciate_emissions(emissions::DataFrame, gspro::DataFrame, gsref::DataFrame;
        basis::Symbol=:mass) -> DataFrame

Apply chemical speciation to convert inventory pollutants into model species.

For each emissions record, looks up the speciation profile via GSREF and applies
the speciation factors from GSPRO. The output has `:species` instead of `:POLID`,
with `ANN_VALUE` split according to speciation factors.

# Arguments
- `emissions::DataFrame`: Must have columns `:FIPS`, `:SCC`, `:POLID`, `:ANN_VALUE`.
  May also have `:COUNTRY`, `:LONGITUDE`, `:LATITUDE`, `:Surrogate`.
- `gspro::DataFrame`: Speciation profiles from [`read_gspro`](@ref).
- `gsref::DataFrame`: Cross-reference from [`read_gsref`](@ref).
- `basis::Symbol=:mass`: `:mass` uses `mass_fraction`, `:mole` uses `split_factor/divisor`.

# Returns
A `DataFrame` with the same columns as input but `:POLID` replaced by `:species`,
and `:ANN_VALUE` scaled by the speciation factor. Each input row may produce
multiple output rows (one per model species).

Records with no matching speciation profile are passed through with
`:species` set to the original `:POLID` and `:ANN_VALUE` unchanged.
"""
function speciate_emissions(
        emissions::DataFrame,
        gspro::DataFrame,
        gsref::DataFrame;
        basis::Symbol = :mass,
    )
    # Identify extra columns to carry through
    extra_cols = Symbol[]
    for col in [:COUNTRY, :LONGITUDE, :LATITUDE, :Surrogate]
        if hasproperty(emissions, col)
            push!(extra_cols, col)
        end
    end

    # Pre-allocate result column vectors
    out_fips = String[]
    out_scc = String[]
    out_species = String[]
    out_ann = Float64[]
    out_extra = Dict{Symbol, Vector{Any}}(col => Any[] for col in extra_cols)

    for erow in eachrow(emissions)
        fips = string(erow.FIPS)
        scc_raw = string(erow.SCC)
        scc = lpad(scc_raw, 10, '0')
        polid = uppercase(string(erow.POLID))
        ann_value = Float64(erow.ANN_VALUE isa Number ? erow.ANN_VALUE : ustrip(erow.ANN_VALUE))

        profile_code = _match_speciation_profile(gsref, fips, scc, polid)

        if profile_code === nothing
            # Pass through unmatched records
            push!(out_fips, fips)
            push!(out_scc, scc_raw)
            push!(out_species, polid)
            push!(out_ann, ann_value)
            for col in extra_cols
                push!(out_extra[col], erow[col])
            end
            continue
        end

        # Find matching entries in GSPRO
        speciated = false
        for prow in eachrow(gspro)
            if prow.profile_code == profile_code && prow.pollutant_id == polid
                if basis == :mole
                    factor = prow.divisor != 0.0 ? prow.split_factor / prow.divisor : 0.0
                else
                    factor = prow.mass_fraction
                end

                push!(out_fips, fips)
                push!(out_scc, scc_raw)
                push!(out_species, prow.species_id)
                push!(out_ann, ann_value * factor)
                for col in extra_cols
                    push!(out_extra[col], erow[col])
                end
                speciated = true
            end
        end

        if !speciated
            # Profile code found but no entries in GSPRO - pass through
            push!(out_fips, fips)
            push!(out_scc, scc_raw)
            push!(out_species, polid)
            push!(out_ann, ann_value)
            for col in extra_cols
                push!(out_extra[col], erow[col])
            end
        end
    end

    result = DataFrame(
        FIPS = out_fips,
        SCC = out_scc,
        species = out_species,
        ANN_VALUE = out_ann,
    )
    for col in extra_cols
        result[!, col] = out_extra[col]
    end

    return result
end
