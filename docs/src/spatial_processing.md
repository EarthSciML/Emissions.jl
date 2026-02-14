# Spatial Processing Pipeline

## Overview

The Emissions.jl package provides a set of mid-level pipeline functions that bridge the gap
between low-level building blocks (FF10 readers, surrogate matrices, spatial indexing) and
the complete spatial processing workflow demonstrated in the
[reference notebook](https://github.com/EarthSciML/Emissions.jl/blob/main/examples/2019neiEmis_3.ipynb).

These pipeline functions handle the **data preparation phase** of emissions processing:
1. Reading multiple FF10-format emissions files
2. Aggregating and filtering emissions records
3. Mapping pollutant identifiers to standard group names
4. Assigning spatial surrogates via grid reference joins with fallback logic
5. Identifying which surrogate shapefiles are actually needed

These functions operate on DataFrames and do not require specific file formats, so users can
load data however they want and pass DataFrames to the pipeline.

**Note**: This pipeline covers the initial data preparation steps. For complete spatial processing
including matrix generation, surrogate computation, and emissions output, users should combine
these functions with the lower-level spatial processing functions (`generate_data_sparse_matrices`,
`generate_weight_sparse_matrices`, etc.) as shown in the examples.

## Pipeline Functions

### File I/O

```@docs
read_ff10
normalize_country
read_gridref
```

### Emissions Processing

```@docs
aggregate_emissions
filter_known_pollutants
map_pollutant_names!
```

### Surrogate Assignment

```@docs
assign_surrogates
build_data_weight_map
find_surrogate_by_code(::Vector{SurrogateSpec}, ::String, ::Int)
```

## Implementation

### Pipeline Workflow Example

```@example pipeline
using Emissions
using DataFrames

# Create synthetic emissions data (simulating two sectors)
nonpoint_data = DataFrame(
    POLID = ["NOX", "VOC", "NOX", "SO2"],
    COUNTRY = ["USA", "USA", "USA", "USA"],
    FIPS = ["36001", "36001", "36005", "36005"],
    SCC = ["2103007000", "2103007000", "2103007000", "2103007000"],
    ANN_VALUE = [150.5, 75.2, 200.1, 125.8]
)

onroad_data = DataFrame(
    POLID = ["NOX", "EXH__VOC", "PM25-PRI"],
    COUNTRY = ["USA", "USA", "USA"],
    FIPS = ["36001", "36001", "36005"],
    SCC = ["2201001000", "2201001000", "2201001000"],
    ANN_VALUE = [50.0, 30.0, 15.0]
)

# Step 1: Aggregate emissions from multiple sectors
combined = aggregate_emissions([nonpoint_data, onroad_data])
println("Combined records: ", nrow(combined))
```

```@example pipeline
# Step 2: Filter to known pollutants only
filtered = filter_known_pollutants(combined)
println("After filtering: ", nrow(filtered))
```

```@example pipeline
# Step 3: Map pollutant names to standard groups
map_pollutant_names!(filtered)
println("Unique pollutants: ", unique(filtered.POLID))
```

```@example pipeline
# Step 4: Assign surrogates
gridref = DataFrame(
    COUNTRY = ["USA", "USA", "USA"],
    FIPS = ["36001", "36005", "00000"],
    SCC = ["2103007000", "2103007000", "2201001000"],
    Surrogate = [100, 100, 200]
)

with_surrogates = assign_surrogates(filtered, gridref)
println("Surrogates assigned: ", count(!ismissing, with_surrogates.Surrogate),
    " of ", nrow(with_surrogates))
```

```@example pipeline
# Step 5: Build the data/weight shapefile map
srgSpecs = [
    SurrogateSpec("USA", "Population", 100, "/data/pop.shp", "POP",
        "/weight/pop_w.shp", "Population surrogate",
        String[], String[], Float64[], "", String[], Float64[]),
    SurrogateSpec("USA", "Roads", 200, "/data/roads.shp", "ROADS",
        "/weight/roads_w.shp", "Road network surrogate",
        String[], String[], Float64[], "", String[], Float64[]),
]

shapefile_map = build_data_weight_map(with_surrogates, srgSpecs)
for (key, labels) in shapefile_map
    println("$(key[1]) + $(key[2]) => $(labels)")
end
```

### Complete Workflow Integration

Here's how these pipeline functions connect to the complete spatial processing workflow:

```@example pipeline
# After using the pipeline functions above, continue with spatial processing:

# Step 6: Create spatial processor configuration
# This requires grid definitions, which would typically be loaded from files
println("Pipeline functions prepare data for spatial processing...")
println("Next steps would involve:")
println("1. Loading grid definitions with read_grid() or NewGridIrregular()")
println("2. Setting up spatial processor with NewSpatialProcessor()")
println("3. Generating sparse matrices with generate_data_sparse_matrices()")
println("4. Computing surrogates with generate_countySurrogate()")
println("5. Writing output with writeEmis()")
```

## Analysis

### Surrogate Assignment with Fallback

The `assign_surrogates` function implements a two-pass matching strategy:
1. **Direct match**: Joins emissions with grid reference on `(COUNTRY, FIPS, SCC)`
2. **Fallback match**: For unmatched records, retries with `FIPS="00000"` (state-level default)

This ensures that emissions records always get a surrogate assignment when a default is available,
even if the specific county-level entry is missing from the grid reference file.

```@example pipeline
# Demonstrate fallback: FIPS 36999 is not in gridref, but "00000" provides a default
test_emissions = DataFrame(
    POLID = ["NOX", "NOX"],
    COUNTRY = ["USA", "USA"],
    FIPS = ["36001", "36999"],
    SCC = ["2103007000", "2103007000"],
    ANN_VALUE = [100.0, 50.0]
)
test_gridref = DataFrame(
    COUNTRY = ["USA", "USA"],
    FIPS = ["36001", "00000"],
    SCC = ["2103007000", "2103007000"],
    Surrogate = [100, 300]
)
result = assign_surrogates(test_emissions, test_gridref)
println("FIPS 36001 surrogate: ", result[result.FIPS .== "36001", :Surrogate][1])
println("FIPS 36999 surrogate (fallback): ", result[result.FIPS .== "36999", :Surrogate][1])
println("FIPS 36999 preserved: ", result[result.FIPS .== "36999", :FIPS][1])
```

### Country Code Normalization

The `normalize_country` function handles the various country code formats found in different
data sources:

```@example pipeline
codes = ["US", "0", "1", "2", "USA", "Canada"]
for code in codes
    println("\"$code\" => \"$(normalize_country(code))\"")
end
```

## Relationship to Reference Workflow

These pipeline functions implement the **data preparation steps 1-6** from the complete
emissions spatial processing workflow demonstrated in the
[2019 NEI emissions notebook](https://github.com/EarthSciML/Emissions.jl/blob/main/examples/2019neiEmis_3.ipynb).

The complete workflow includes additional steps:
- **Steps 7-9**: Sparse matrix generation from shapefiles (`generate_data_sparse_matrices`, `generate_weight_sparse_matrices`, `generate_grid_sparse_matrices`)
- **Steps 10-11**: Surrogate computation and location index refinement (`generate_countySurrogate`, `update_locIndex`)
- **Step 12**: Emissions output to shapefiles (`writeEmis`)

These additional functions are available in Emissions.jl but operate at a lower level,
requiring more detailed configuration of grid definitions and spatial processors.
