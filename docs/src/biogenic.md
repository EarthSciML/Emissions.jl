# Biogenic Emissions

## Overview

Biogenic emissions (isoprene, terpenes, etc.) from vegetation are computed using
the BEIS (Biogenic Emission Inventory System) methodology, implementing
the SMOKE Normbeis functionality.

The computation combines:
- Land cover data (BELD format)
- Vegetation-specific emission factors
- Temperature and light response functions (Guenther, 1993)

```@docs
BiogenicConfig
read_beld
read_emission_factors
temperature_adjustment
light_adjustment
compute_biogenic_emissions
```

## Usage

### Basic Computation

```julia
using Emissions

# Configure
config = BiogenicConfig(
    "path/to/beld.csv",
    "path/to/emission_factors.csv",
    :summer
)

# Define grid and meteorology
grid = NewGridIrregular("test", 10, 10, "EPSG:4326", 0.5, 0.5, -80.0, 35.0)
temperature = fill(303.15, 100)  # K per grid cell
par = fill(1000.0, 100)          # μmol/m²/s per grid cell

# Compute
biogenic = compute_biogenic_emissions(config, grid, temperature, par)
```

### Combining with Anthropogenic Emissions

Biogenic output is compatible with `merge_categories`:

```julia
# anthropogenic = process_emissions(...)  # gridded hourly
# biogenic = compute_biogenic_emissions(...)
# combined = merge_categories(anthropogenic, biogenic)
```

## Temperature and Light Response

### Isoprene

- **Temperature**: Arrhenius-type function with optimum near 314 K
- **Light**: Hyperbolic PAR response (Guenther, 1993)
- Both temperature and PAR adjustments are applied

### Terpenes and Other BVOCs

- **Temperature**: Exponential response (β = 0.09 K⁻¹)
- **Light**: No light dependence (emissions are temperature-driven only)

## File Formats

### BELD Land Cover File

Comma-delimited: `cell_index,land_use_type,fraction`

### Emission Factors File

Comma-delimited: `land_use_type,species,summer_factor,winter_factor`

Factors are in μg/m²/hr at standard conditions (30°C, 1000 μmol/m²/s PAR).
