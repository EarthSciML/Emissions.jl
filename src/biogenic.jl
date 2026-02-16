export BiogenicConfig, read_beld, read_emission_factors,
    temperature_adjustment, light_adjustment, compute_biogenic_emissions

"""
    BiogenicConfig

Configuration for biogenic emissions computation.

# Fields
- `land_use_file::String`: Path to BELD land cover data file
- `emission_factors_file::String`: Path to emission factor data file
- `season::Symbol`: Season (`:summer` or `:winter`) for factor selection
"""
struct BiogenicConfig
    land_use_file::String
    emission_factors_file::String
    season::Symbol
end

"""
    read_beld(filepath::AbstractString, grid::GridDef) -> DataFrame

Read BELD (Biogenic Emissions Landcover Database) land cover data.

Lines starting with `#` are skipped. Each data line has comma-delimited fields:
`cell_index,land_use_type,fraction`

# Returns
A `DataFrame` with columns:
- `:cell_index` (Int): Grid cell index (1-based)
- `:land_use_type` (String): Land use category name
- `:fraction` (Float64): Fraction of cell covered by this land use (0-1)
"""
function read_beld(filepath::AbstractString, grid::GridDef)
    lines = readlines(filepath)
    records = DataFrame(
        cell_index = Int[],
        land_use_type = String[],
        fraction = Float64[],
    )
    ncells = grid.Nx * grid.Ny
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        # Remove inline comments
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ',')
        if length(parts) >= 3
            cell_index = parse(Int, strip(parts[1]))
            land_use_type = strip(parts[2])
            fraction = parse(Float64, strip(parts[3]))
            if 1 <= cell_index <= ncells && fraction > 0.0
                push!(
                    records, (
                        cell_index = cell_index,
                        land_use_type = land_use_type,
                        fraction = fraction,
                    )
                )
            end
        end
    end
    return records
end

"""
    read_emission_factors(filepath::AbstractString) -> DataFrame

Read biogenic emission factor data (B3FAC/B4FAC format).

Lines starting with `#` are skipped. Each data line has comma-delimited fields:
`land_use_type,species,summer_factor,winter_factor`

Emission factors are in units of μg/m²/hr at standard conditions (30°C, 1000 μmol/m²/s PAR).

# Returns
A `DataFrame` with columns:
- `:land_use_type` (String): Land use category name
- `:species` (String): Biogenic species name (e.g., "ISOP", "TERP", "NO")
- `:summer_factor` (Float64): Summer emission factor (μg/m²/hr)
- `:winter_factor` (Float64): Winter emission factor (μg/m²/hr)
"""
function read_emission_factors(filepath::AbstractString)
    lines = readlines(filepath)
    records = DataFrame(
        land_use_type = String[],
        species = String[],
        summer_factor = Float64[],
        winter_factor = Float64[],
    )
    for line in lines
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        if occursin('!', line)
            line = strip(split(line, '!')[1])
        end
        parts = split(line, ',')
        if length(parts) >= 4
            land_use_type = strip(parts[1])
            species = uppercase(strip(parts[2]))
            summer_factor = parse(Float64, strip(parts[3]))
            winter_factor = parse(Float64, strip(parts[4]))
            push!(
                records, (
                    land_use_type = land_use_type,
                    species = species,
                    summer_factor = summer_factor,
                    winter_factor = winter_factor,
                )
            )
        end
    end
    return records
end

"""
    temperature_adjustment(temp_K::Real, species::AbstractString) -> Float64

Compute the Guenther (1993) temperature response factor for biogenic emissions.

For isoprene (`"ISOP"`), uses the Arrhenius-type function:
`γ_T = CT1 * exp(CT2 * (T - Ts) / (R * T * Ts)) / (1 + exp(CT3 * (T - TM) / (R * T * Ts)))`

For terpenes and other species, uses an exponential function:
`γ_T = exp(β * (T - Ts))`

where `Ts = 303.15 K` is the standard temperature.

# Arguments
- `temp_K`: Temperature in Kelvin
- `species`: Species name ("ISOP" for isoprene, anything else for terpene-like response)

# Returns
Temperature adjustment factor (dimensionless, 1.0 at standard conditions for terpenes).
"""
function temperature_adjustment(temp_K::Real, species::AbstractString)
    T = Float64(temp_K)
    Ts = 303.15  # Standard temperature (30°C in K)
    R = 8.314    # Universal gas constant (J/mol/K)

    if uppercase(species) == "ISOP"
        # Guenther (1993) isoprene temperature response
        CT1 = 95000.0   # J/mol
        CT2 = 230000.0  # J/mol
        TM = 314.0      # K (optimum temperature)

        numerator = exp(CT1 * (T - Ts) / (R * T * Ts))
        denominator = 1.0 + exp(CT2 * (T - TM) / (R * T * Ts))
        return max(numerator / denominator, 0.0)
    else
        # Exponential response for terpenes and other BVOCs
        β = 0.09  # Temperature sensitivity (1/K)
        return exp(β * (T - Ts))
    end
end

"""
    light_adjustment(par::Real) -> Float64

Compute the Guenther (1993) PAR response factor for isoprene emissions.

`γ_L = α * CL1 * PAR / sqrt(1 + α² * PAR²)`

where PAR is in μmol/m²/s, `α = 0.0027`, `CL1 = 1.066`.

# Arguments
- `par`: Photosynthetically Active Radiation in μmol/m²/s

# Returns
Light adjustment factor (dimensionless, 1.0 at PAR = 1000 μmol/m²/s).
"""
function light_adjustment(par::Real)
    PAR = Float64(par)
    if PAR <= 0.0
        return 0.0
    end
    α = 0.0027
    CL1 = 1.066
    return α * CL1 * PAR / sqrt(1.0 + α^2 * PAR^2)
end

"""
    compute_biogenic_emissions(config::BiogenicConfig, grid::GridDef,
        temperature::Vector{Float64}, par::Vector{Float64}) -> DataFrame

Compute gridded biogenic emissions for a single timestep.

Combines land use data, emission factors, and meteorological adjustments to
produce gridded biogenic emissions compatible with [`merge_categories`](@ref).

# Arguments
- `config::BiogenicConfig`: Configuration with paths and season.
- `grid::GridDef`: Target grid definition.
- `temperature::Vector{Float64}`: Temperature per grid cell (K). Length must equal `grid.Nx * grid.Ny`.
- `par::Vector{Float64}`: PAR per grid cell (μmol/m²/s). Length must equal `grid.Nx * grid.Ny`.

# Returns
A `DataFrame` with columns:
- `:grid_row` (Int): Grid row index
- `:grid_col` (Int): Grid column index
- `:species` (String): Biogenic species name
- `:emission_rate` (Float64): Emission rate (μg/m²/hr adjusted)
"""
function compute_biogenic_emissions(
        config::BiogenicConfig,
        grid::GridDef,
        temperature::Vector{Float64},
        par::Vector{Float64},
    )
    ncells = grid.Nx * grid.Ny
    length(temperature) == ncells || throw(
        ArgumentError("temperature vector length ($(length(temperature))) must match grid cells ($ncells)")
    )
    length(par) == ncells || throw(
        ArgumentError("par vector length ($(length(par))) must match grid cells ($ncells)")
    )

    # Read input data
    land_use = read_beld(config.land_use_file, grid)
    ef = read_emission_factors(config.emission_factors_file)

    # Select season factor column
    factor_col = config.season == :winter ? :winter_factor : :summer_factor

    # Build lookup: (land_use_type, species) -> emission factor
    ef_lookup = Dict{Tuple{String, String}, Float64}()
    all_species = unique(ef.species)
    for row in eachrow(ef)
        ef_lookup[(row.land_use_type, row.species)] = row[factor_col]
    end

    # Accumulate emissions: (cell_index, species) -> rate
    accumulator = Dict{Tuple{Int, String}, Float64}()

    for lu_row in eachrow(land_use)
        cell_idx = lu_row.cell_index
        lu_type = lu_row.land_use_type
        frac = lu_row.fraction
        T = temperature[cell_idx]
        P = par[cell_idx]

        for sp in all_species
            base_ef = get(ef_lookup, (lu_type, sp), 0.0)
            if base_ef <= 0.0
                continue
            end

            # Apply temperature and light adjustments
            γ_T = temperature_adjustment(T, sp)
            if uppercase(sp) == "ISOP"
                γ_L = light_adjustment(P)
                rate = base_ef * frac * γ_T * γ_L
            else
                rate = base_ef * frac * γ_T
            end

            if rate > 0.0
                acc_key = (cell_idx, sp)
                accumulator[acc_key] = get(accumulator, acc_key, 0.0) + rate
            end
        end
    end

    # Convert to DataFrame
    n = length(accumulator)
    grid_rows = Vector{Int}(undef, n)
    grid_cols = Vector{Int}(undef, n)
    species_out = Vector{String}(undef, n)
    rates = Vector{Float64}(undef, n)

    for (idx, ((cell_idx, sp), rate)) in enumerate(accumulator)
        j = (cell_idx - 1) ÷ grid.Nx + 1
        i = (cell_idx - 1) % grid.Nx + 1
        grid_rows[idx] = j
        grid_cols[idx] = i
        species_out[idx] = sp
        rates[idx] = rate
    end

    result = DataFrame(
        grid_row = grid_rows,
        grid_col = grid_cols,
        species = species_out,
        emission_rate = rates,
    )
    sort!(result, [:grid_row, :grid_col, :species])
    return result
end
