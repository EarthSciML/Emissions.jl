# Emissions Controls

## Overview

Emissions controls implement the SMOKE Cntlmat functionality for applying growth
factors and control efficiencies to scale emissions. This enables scenario analysis,
future-year projections, and regulatory compliance modeling.

```@docs
ControlSpec
read_growth_factors
read_control_factors
apply_controls
```

## Usage

### Growth Factors

Growth factors project emissions from a base year to a target year:

```julia
using Emissions

controls = read_growth_factors("path/to/growth_factors.txt")
projected = apply_controls(emissions, controls)
```

### Control Factors

Control factors apply efficiency/effectiveness/penetration reductions:

```julia
controls = read_control_factors("path/to/control_factors.txt")
controlled = apply_controls(emissions, controls)
```

### Combined Controls

Growth and multiplicative controls can be combined. Growth is applied first,
then multiplicative:

```julia
growth = read_growth_factors("growth.txt")
mult = read_control_factors("control.txt")
all_controls = vcat(growth, mult)
result = apply_controls(emissions, all_controls)
```

## File Formats

### Growth Factor File

Semicolon-delimited: `FIPS;SCC;pollutant;base_year;target_year;growth_factor`

### Control Factor File

Semicolon-delimited: `FIPS;SCC;pollutant;efficiency;effectiveness;penetration`

The multiplicative factor is: `1 - (efficiency/100 * effectiveness/100 * penetration/100)`

## Hierarchical Matching

Control matching uses a 7-level hierarchy from most specific to least:
1. Exact FIPS + SCC + pollutant
2. Exact FIPS + SCC (any pollutant)
3. Exact FIPS + pollutant (any SCC)
4. Exact FIPS only
5. National default + SCC + pollutant
6. National default + pollutant
7. National default only
