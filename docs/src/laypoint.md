# Vertical Layer Allocation

## Overview

Vertical layer allocation implements the SMOKE Laypoint functionality. It distributes
elevated point source emissions across vertical atmospheric layers using meteorological
data and the ASME (1973) plume rise calculation.

```@docs
LayerConfig
MetProfile
compute_layer_fractions
allocate_point_to_layers
laypoint
```

## Usage

### Basic Layer Allocation

```julia
using Emissions

# Define atmospheric layers (staggered heights in meters)
config = LayerConfig([0.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0])

# Define meteorological profile
met = MetProfile(
    fill(280.0, 5),    # temperature per layer (K)
    fill(5.0, 5),      # wind speed per layer (m/s)
    fill(0.0, 5),      # stability class (0=unstable)
    fill(0.0, 5),      # stability parameter
)
met_profiles = Dict("default" => met)

# Allocate point sources to layers
result = laypoint(point_emissions, met_profiles, config)
```

### Integration with Merge

The output from `laypoint` includes `:layer` and `:layer_fraction` columns
that are automatically handled by `merge_emissions` when present, producing
3D gridded emissions with layer information.

## Algorithm

1. For each point source, compute ASME (1973) plume rise using stack
   parameters and meteorological data
2. Determine plume bottom (stack height) and plume top (stack + rise + spread)
3. Distribute the plume across layers proportionally to overlap height
4. Layer fractions sum to 1.0 for each source
5. Plumes above model top are allocated to the top layer
