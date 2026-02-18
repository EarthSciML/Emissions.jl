# QA Reporting

## Overview

The QA reporting module provides summary statistics and comparison tools for
emissions data, analogous to SMOKE's Smkreport utility. Reports can be generated
from both pre-merge inventory data and post-merge gridded emissions.

```@docs
ReportConfig
summary_by_pollutant
summary_by_region
summary_by_scc
summary_by_time
compare_inventories
emissions_report
```

## Usage

### Summary by Pollutant

```@example reporting
using Emissions, DataFrames

# Pre-merge inventory data
inventory = DataFrame(
    FIPS = ["36001", "36001", "36005", "06001"],
    SCC = ["2103007000", "2103007000", "2103007001", "2103007000"],
    POLID = ["NOX", "NOX", "VOC", "NOX"],
    ANN_VALUE = [100.0, 50.0, 200.0, 75.0]
)

summary_by_pollutant(inventory)
```

### Summary by SCC

```@example reporting
summary_by_scc(inventory)
```

### Inventory Comparison

```@example reporting
# Compare base and scenario inventories
base = DataFrame(
    FIPS = ["36001", "36005"],
    SCC = ["2103007000", "2103007001"],
    POLID = ["NOX", "VOC"],
    ANN_VALUE = [100.0, 200.0]
)

scenario = DataFrame(
    FIPS = ["36001", "36005"],
    SCC = ["2103007000", "2103007001"],
    POLID = ["NOX", "VOC"],
    ANN_VALUE = [80.0, 210.0]
)

compare_inventories(base, scenario)
```

### Custom Reports

```@example reporting
config = ReportConfig(group_by = [:FIPS, :POLID], top_n = 5)
emissions_report(inventory; config = config)
```

### Time Series Reports (Post-Merge)

```@example reporting
using Dates

gridded = DataFrame(
    grid_row = [1, 1, 1, 1],
    grid_col = [1, 1, 1, 1],
    hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 12),
            DateTime(2019, 7, 2, 0), DateTime(2019, 7, 2, 12)],
    pollutant = ["NOX", "NOX", "NOX", "NOX"],
    emission_rate = [100.0, 150.0, 90.0, 140.0]
)

summary_by_time(gridded; resolution = :daily)
```
