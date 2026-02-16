# Elevated Source Identification

## Overview

Elevated source identification implements the SMOKE Elevpoint functionality.
It classifies point sources as surface-level, elevated, or Plume-in-Grid (PinG)
based on stack parameters and analytical plume rise estimates. This classification
determines how sources are treated in vertical layer allocation and air quality
modeling.

```@docs
ElevationCriteria
analytical_plume_rise
classify_point_sources
group_stacks
```

## Usage

### Classifying Point Sources

```julia
using Emissions

# Classify with default SMOKE criteria
classified = classify_point_sources(point_emissions)

# Filter by class
elevated = filter(r -> r.source_class == "elevated", classified)
ping = filter(r -> r.source_class == "ping", classified)
surface = filter(r -> r.source_class == "surface", classified)
```

### Custom Criteria

```julia
criteria = ElevationCriteria(
    20.0,   # min stack height (m)
    5.0,    # min exit velocity (m/s)
    400.0,  # min exit temperature (K)
    5.0,    # min flow rate (m³/s)
    50.0,   # min plume rise (m)
    50.0,   # PinG stack height (m)
    0.01,   # PinG emissions threshold (kg/s)
)
classified = classify_point_sources(point_emissions; criteria=criteria)
```

### Stack Grouping

```julia
# Group similar stacks at same facility
grouped = group_stacks(point_emissions; height_bin=10.0, temp_bin=50.0)
```

## Classification Logic

A source is **elevated** if ANY of these criteria are met:
- Stack height ≥ threshold
- Exit velocity ≥ threshold
- Exit temperature ≥ threshold
- Flow rate ≥ threshold
- Analytical plume rise ≥ threshold

An elevated source is further classified as **PinG** if ALL PinG criteria are met:
- Stack height ≥ PinG stack height threshold
- Emissions rate ≥ PinG emissions threshold (if set)
