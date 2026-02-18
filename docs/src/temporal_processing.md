# Temporal Processing

## Overview

The temporal allocation module converts annual emissions inventory data to hourly emissions
using temporal profiles. This implements the SMOKE `Temporal` program functionality.

Temporal allocation uses three types of profiles:
- **Monthly profiles**: 12 factors distributing annual total across months (sum to 1.0)
- **Weekly profiles**: 7 factors for day-of-week weighting (sum to 7.0)
- **Diurnal profiles**: 24 factors distributing daily total across hours (sum to 1.0)

## API Reference

### Profile I/O

```@docs
read_temporal_profiles
read_temporal_xref
read_day_specific
read_hour_specific
```

### Temporal Allocation

```@docs
temporal_allocate
```

### Merging

After temporal and spatial allocation, emissions from different source categories
are merged into final gridded hourly output.

```@docs
merge_emissions
merge_categories
```

## Example

```@example temporal
using Emissions
using DataFrames
using Dates

# Create synthetic annual emissions
emissions = DataFrame(
    FIPS = ["36001", "36005"],
    SCC = ["2103007000", "2103007000"],
    POLID = ["NOX", "VOC"],
    ANN_VALUE = [100.0, 50.0]
)

# Create uniform temporal profiles
profiles = DataFrame(
    profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
    profile_id = [1, 1, 1],
    factors = [fill(1.0/12.0, 12), fill(1.0, 7), fill(1.0/24.0, 24)]
)

# Create cross-reference mapping all SCCs to profile ID 1
xref = DataFrame(
    FIPS = ["00000"],
    SCC = ["2103007000"],
    monthly_id = [1],
    weekly_id = [1],
    diurnal_id = [1]
)

# Allocate for a 24-hour episode
ep_start = DateTime(2019, 7, 1, 0)
ep_end = DateTime(2019, 7, 2, 0)

hourly = temporal_allocate(emissions, profiles, xref, ep_start, ep_end)
println("Hourly records: ", nrow(hourly))
first(hourly, 5)
```

```@example temporal
# Merge with spatial allocation
using SparseArrays

grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
locIndex = Dict{String, IndexInfo}(
    "36001" => IndexInfo([1, 1], [1, 2], [0.6, 0.4], true, true),
    "36005" => IndexInfo([2], [1], [1.0], true, true)
)

gridded = merge_emissions(hourly, locIndex, grid)
println("Gridded records: ", nrow(gridded))
first(gridded, 5)
```
