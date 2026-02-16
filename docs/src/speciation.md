# Chemical Speciation

## Overview

Chemical speciation converts inventory-level pollutants (e.g., VOC, NOX) into
model-ready chemical species (e.g., FORM, ALD2, NO, NO2) required by air quality
models such as CMAQ and CAMx. This implements the SMOKE Spcmat functionality.

Speciation uses two input files:
- **GSPRO**: Speciation profile file containing split factors and mass fractions
- **GSREF**: Cross-reference file mapping sources to speciation profiles

```@docs
read_gspro
read_gsref
build_speciation_matrix
speciate_emissions
```

## Usage

### Reading Speciation Files

```julia
using Emissions

# Read speciation profiles and cross-reference
gspro = read_gspro("path/to/gspro.txt")
gsref = read_gsref("path/to/gsref.txt")
```

### Applying Speciation

```julia
# Apply mass-based speciation to emissions
speciated = speciate_emissions(emissions, gspro, gsref; basis=:mass)

# Or mole-based speciation
speciated_mole = speciate_emissions(emissions, gspro, gsref; basis=:mole)
```

### Pipeline Integration

Speciation integrates into the emissions processing pipeline after pollutant
name mapping and before spatial allocation:

```julia
result = process_emissions(
    inventory_files = [("inventory.csv", :nonpoint)],
    grid = grid,
    gspro = gspro,
    gsref = gsref,
    speciation_basis = :mass,
    # ... other arguments
)
```

## File Formats

### GSPRO Format

Semicolon-delimited with 6 fields per line:
```
profile_code;pollutant_id;species_id;split_factor;divisor;mass_fraction
```

### GSREF Format

Semicolon-delimited with 4 fields per line:
```
FIPS;SCC;pollutant_id;profile_code
```

### Hierarchical Matching

GSREF matching follows a 3-level hierarchy:
1. Exact FIPS + SCC + pollutant match
2. National default (FIPS="00000") + SCC + pollutant
3. Pollutant-only match (SCC="0000000000")
