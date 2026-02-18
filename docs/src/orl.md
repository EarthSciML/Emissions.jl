# ORL Format

## Overview

The ORL (One Record per Line) format is a legacy EPA inventory format used by
SMOKE prior to the FF10 format. This module provides readers for all ORL
inventory types, following the same pattern as the FF10 readers.

Supported ORL formats:
- **Nonpoint** (ARINV): Area source emissions
- **Point** (PTINV): Point source emissions with stack parameters
- **Nonroad** (ARINV): Mobile nonroad source emissions
- **Onroad** (MBINV): Mobile on-road source emissions
- **Fire** (PTFIRE): Wildfire and prescribed fire emissions

```@docs
ORLNonPointDataFrame
ORLPointDataFrame
ORLNonRoadDataFrame
ORLOnRoadDataFrame
ORLFireDataFrame
read_orl
```

## Usage

### Reading ORL Files

```@example orl
using Emissions

# Read different ORL format types
# emis = read_orl("nonpoint_inventory.csv", :nonpoint)
# emis = read_orl("point_inventory.csv", :point)
# emis = read_orl("fire_inventory.csv", :fire)
```

### Unit Conversions

All ORL readers automatically perform the following conversions:
- **Emission values**: tons/year to kg/s (using `tonperyear` conversion factor)
- **FIPS codes**: Standardized to 5-digit format
- **SCC codes**: Left-padded to 10 digits

For point sources, stack parameters are additionally converted:
- Stack height: feet to meters
- Stack diameter: feet to meters
- Stack temperature: Fahrenheit to Kelvin
- Stack flow: ft³/s to m³/s
- Stack velocity: ft/s to m/s

### Pipeline Integration

ORL files can be used in the processing pipeline by setting `inventory_type=:orl`:

```julia
result = process_emissions(;
    inventory_files = [("orl_inventory.csv", :nonpoint)],
    grid = my_grid,
    inventory_type = :orl,  # Use ORL format readers
    # ... other options
)
```
