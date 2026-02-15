# Spatial Processing Pipeline

## Overview

The Emissions.jl package provides a complete spatial processing pipeline for allocating
emissions to grid cells, implementing the workflow from the
[2019 NEI emissions processing notebook](https://github.com/EarthSciML/Emissions.jl/blob/main/examples/2019neiEmis_3.ipynb).

The pipeline consists of two layers:

1. **Data Preparation**: Reading FF10 files, aggregating emissions, filtering pollutants,
   mapping names, and assigning spatial surrogates.
2. **Spatial Allocation**: Computing grid cell indices for emission locations, optionally
   refining with surrogate-weighted fractions, and distributing emissions to grid cells.

## Data Preparation Functions

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

## Spatial Allocation Functions

### High-Level Workflow Functions

```@docs
location_key
compute_grid_indices
refine_indices_with_surrogates
allocate_emissions_to_grid
process_emissions_spatial
```

### Low-Level Surrogate Generation Functions

The following low-level functions are used for advanced surrogate generation from shapefiles.
For complete API documentation, see the [NEI Processing](nei_processing.md#surrogate-operations) page.

- [`generate_data_sparse_matrices`](@ref) - Generate sparse matrices from data shapefiles
- [`generate_weight_sparse_matrices`](@ref) - Generate weight matrices from surrogate shapefiles
- [`generate_grid_sparse_matrices`](@ref) - Generate grid area matrices
- [`generate_countySurrogate`](@ref) - Combine data and weight matrices into normalized surrogates
- [`update_locIndex`](@ref) - Update location indices with surrogate data

## Implementation

### Complete Workflow Example

The following example demonstrates the full spatial processing pipeline using
synthetic data. This corresponds to the workflow in the reference notebook.

```@example workflow
using Emissions
using DataFrames
using SparseArrays

# =============================================================================
# Step 1: Create Synthetic Emissions Data
# =============================================================================
# Simulating two emission sectors with point source coordinates
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

println("Nonpoint records: ", nrow(nonpoint_data))
println("Onroad records: ", nrow(onroad_data))
```

```@example workflow
# =============================================================================
# Step 2: Aggregate Emissions from Multiple Sectors
# =============================================================================
combined = aggregate_emissions([nonpoint_data, onroad_data])
println("Combined records: ", nrow(combined))
first(combined, 5)
```

```@example workflow
# =============================================================================
# Step 3: Filter to Known Pollutants and Map Names
# =============================================================================
filtered = filter_known_pollutants(combined)
println("After filtering: ", nrow(filtered))

map_pollutant_names!(filtered)
println("Unique pollutants: ", unique(filtered.POLID))
```

```@example workflow
# =============================================================================
# Step 4: Assign Spatial Surrogates
# =============================================================================
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

```@example workflow
# =============================================================================
# Step 5: Create Target Grid
# =============================================================================
# Create a simple 4x4 grid (in real use, this would be an InMAP or CMAQ grid)
grid = NewGridIrregular("CONUS_4x4", 4, 4, "EPSG:4326", 1.0, 1.0, -76.0, 39.0)
println("Grid: $(grid.Nx) x $(grid.Ny) = $(grid.Nx * grid.Ny) cells")
```

```@example workflow
# =============================================================================
# Step 6: Spatial Allocation with Surrogates
# =============================================================================
# Create synthetic county surrogates that distribute emissions across grid cells.
# In real use, these are generated by generate_data_sparse_matrices(),
# generate_weight_sparse_matrices(), and generate_countySurrogate().

# County 36001 (Albany, NY): 60% in cell (1,1), 40% in cell (1,2)
srg_36001 = sparse([1, 1], [1, 2], [0.6, 0.4], 4, 4)

# County 36005 (Bronx, NY): 50% in cell (2,1), 30% in cell (2,2), 20% in cell (3,1)
srg_36005 = sparse([2, 2, 3], [1, 2, 1], [0.5, 0.3, 0.2], 4, 4)

county_surrogates = Dict("36001" => srg_36001, "36005" => srg_36005)

# Run the complete spatial allocation workflow
gridded = process_emissions_spatial(with_surrogates, grid;
    county_surrogates=county_surrogates)

println("Gridded output: $(nrow(gridded)) rows")
gridded
```

### Step-by-Step Spatial Allocation

The `process_emissions_spatial` function orchestrates three steps that can
also be called individually for more control:

```@example workflow
# Step 6a: Compute grid indices for all unique locations
locIndex = compute_grid_indices(with_surrogates, grid)
println("Location indices computed: ", length(locIndex))
for (key, idx) in locIndex
    println("  $key: inGrid=$(idx.inGrid), cells=$(length(idx.rows))")
end
```

```@example workflow
# Step 6b: Refine area-source indices with surrogate data
refine_indices_with_surrogates(locIndex, county_surrogates)
println("\nAfter surrogate refinement:")
for (key, idx) in locIndex
    println("  $key: inGrid=$(idx.inGrid), cells=$(length(idx.rows))")
    if idx.inGrid
        for k in eachindex(idx.rows)
            println("    cell($(idx.rows[k]),$(idx.cols[k])): $(round(idx.fracs[k]*100, digits=1))%")
        end
    end
end
```

```@example workflow
# Step 6c: Allocate emissions to grid cells
gridded_manual = allocate_emissions_to_grid(with_surrogates, locIndex, grid)
println("Manual allocation result: $(nrow(gridded_manual)) rows")
gridded_manual
```

### Point Source Processing

Point sources with coordinates are allocated directly to their containing grid cell:

```@example workflow
# Point source emissions with explicit coordinates
point_emissions = DataFrame(
    FIPS = ["36001", "36001"],
    POLID = ["NOX", "VOC"],
    SCC = ["2103007000", "2103007000"],
    ANN_VALUE = [1.0e-3, 5.0e-4],
    LONGITUDE = [-75.5, -75.5],
    LATITUDE = [39.5, 39.5]
)

point_result = process_emissions_spatial(point_emissions, grid)
println("Point source allocation:")
point_result
```

## Analysis

### Surrogate Assignment with Fallback

The `assign_surrogates` function implements a two-pass matching strategy:
1. **Direct match**: Joins emissions with grid reference on `(COUNTRY, FIPS, SCC)`
2. **Fallback match**: For unmatched records, retries with `FIPS="00000"` (state-level default)

This ensures that emissions records always get a surrogate assignment when a default
is available, even if the specific county-level entry is missing from the grid reference.

```@example workflow
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

### Emissions Conservation

A key property of the spatial allocation is that total emissions are conserved.
The surrogate fractions are normalized to sum to 1.0 for each county, ensuring
that distributing emissions across grid cells preserves the total:

```@example workflow
# Verify conservation: total emissions in == total emissions out
input_nox = sum(filter(r -> r.POLID == "NOX", with_surrogates).ANN_VALUE)
output_nox = sum(gridded.NOX)
println("Input NOX total:  ", round(input_nox, digits=4))
println("Output NOX total: ", round(output_nox, digits=4))
println("Conservation: ", isapprox(input_nox, output_nox, rtol=1e-10) ? "PASS" : "FAIL")
```

### Country Code Normalization

The `normalize_country` function handles the various country code formats found in
different data sources:

```@example workflow
codes = ["US", "0", "1", "2", "USA", "Canada"]
for code in codes
    println("\"$code\" => \"$(normalize_country(code))\"")
end
```

### Advanced Surrogate Generation

For more complex spatial allocation scenarios, you can generate surrogate matrices directly from shapefiles using the low-level surrogate generation functions. This example demonstrates how to create county surrogates from census and employment data:

```@example workflow
# =============================================================================
# Step 7: Advanced Surrogate Generation from Shapefiles
# =============================================================================
# This demonstrates the workflow for generating surrogates when you have
# county polygons and weight data (e.g., population, employment) shapefiles

println("Advanced surrogate generation workflow:")
println("1. Generate county polygon matrices using:")
println("   data_matrices = generate_data_sparse_matrices(")
println("       \"counties.shp\", \"FIPS\", grid, \"EPSG:4326\")")
println()
println("2. Generate weight matrices (e.g., from census blocks) using:")
println("   weight_matrix = generate_weight_sparse_matrices(")
println("       \"census_blocks.shp\", [\"POPULATION\"], [1.0], grid, \"EPSG:4326\")")
println()
println("3. Generate grid matrices for normalization using:")
println("   grid_matrix = generate_grid_sparse_matrices(grid)")
println()
println("4. Combine into county-level surrogates using:")
println("   county_surrogates = generate_countySurrogate(")
println("       data_matrices, weight_matrix, grid_matrix)")
println()

# For demonstration, we'll show the expected data structure:
println("Example surrogate matrix structure:")
example_matrix = sparse([1, 1, 2], [1, 2, 1], [0.6, 0.4, 1.0], 4, 4)
println("Matrix dimensions: ", size(example_matrix))
println("Non-zero entries: ", nnz(example_matrix))
println("Row sums: ", [sum(example_matrix[i, :]) for i in 1:size(example_matrix, 1)])
```

```@example workflow
# =============================================================================
# Step 8: Complete End-to-End Workflow with Manual Surrogate Creation
# =============================================================================
# This demonstrates creating surrogates programmatically without shapefiles
# (useful for testing and when actual shapefile data is not available)

println("Creating synthetic surrogate data for complete workflow demonstration:")

# Simulate data matrices (county polygons)
# In real use, these come from generate_data_sparse_matrices()
println("1. Simulating data matrices (county coverage):")
data_matrices = Dict{String, SparseMatrixCSC{Float64,Int}}()
data_matrices["36001"] = sparse([1, 1, 2], [1, 2, 2], [0.8, 0.3, 0.5], 4, 4)  # Albany county
data_matrices["36005"] = sparse([2, 3, 3], [1, 1, 2], [0.7, 0.4, 0.6], 4, 4)  # Bronx county
println("   County 36001 covers ", nnz(data_matrices["36001"]), " grid cells")
println("   County 36005 covers ", nnz(data_matrices["36005"]), " grid cells")

# Simulate weight matrix (population distribution)
# In real use, this comes from generate_weight_sparse_matrices()
println("\n2. Simulating weight matrix (population distribution):")
weight_matrix = sparse([1, 1, 2, 2, 3, 3], [1, 2, 1, 2, 1, 2],
                      [1000.0, 800.0, 1500.0, 1200.0, 600.0, 900.0], 4, 4)
println("   Population weight matrix has ", nnz(weight_matrix), " populated cells")

# Generate grid matrix
# In real use, this comes from generate_grid_sparse_matrices()
println("\n3. Generating grid matrix (cell areas):")
grid_matrix = generate_grid_sparse_matrices(grid)
println("   Grid matrix covers ", nnz(grid_matrix), " cells")

# Generate county surrogates using the implemented function
println("\n4. Generating county surrogates:")
county_surrogates_full = generate_countySurrogate(data_matrices, weight_matrix, grid_matrix)
println("   Generated surrogates for ", length(county_surrogates_full), " counties")

for (county, surrogate) in county_surrogates_full
    total = sum(surrogate)
    println("   County $county: $(nnz(surrogate)) cells, total weight = $(round(total, digits=6))")
end

# Now use these surrogates in the workflow
println("\n5. Running spatial allocation with generated surrogates:")
final_result = process_emissions_spatial(with_surrogates, grid;
    county_surrogates=county_surrogates_full)

println("Final gridded emissions:")
final_result
```

## Relationship to Reference Workflow

These functions implement the complete emissions spatial processing workflow from the
[2019 NEI emissions notebook](https://github.com/EarthSciML/Emissions.jl/blob/main/examples/2019neiEmis_3.ipynb).

| Notebook Step | Function | Notes |
|:---|:---|:---|
| Read FF10 files | [`read_ff10`](@ref) | ✓ Complete |
| Concatenate and aggregate | [`aggregate_emissions`](@ref) | ✓ Complete |
| Filter pollutants | [`filter_known_pollutants`](@ref) | ✓ Complete |
| Map pollutant names | [`map_pollutant_names!`](@ref) | ✓ Complete |
| Assign surrogates via grid ref | [`assign_surrogates`](@ref) | ✓ Complete |
| Build shapefile map | [`build_data_weight_map`](@ref) | ✓ Complete |
| Compute grid indices | [`compute_grid_indices`](@ref) | ✓ Complete |
| Generate sparse matrices | [`generate_data_sparse_matrices`](@ref), [`generate_weight_sparse_matrices`](@ref), [`generate_grid_sparse_matrices`](@ref) | ✓ Complete (uses GridDef instead of explicit bounds) |
| Compute county surrogates | [`generate_countySurrogate`](@ref) | ✓ Complete |
| Refine with surrogates | [`refine_indices_with_surrogates`](@ref) | ✓ Complete |
| Allocate to grid | [`allocate_emissions_to_grid`](@ref) | ✓ Complete |
| Complete workflow | [`process_emissions_spatial`](@ref) | ✓ Complete |
| Write output | [`writeEmis`](@ref) | ⚠ Simplified API, outputs DataFrame instead of shapefile |

**Notes on Implementation Differences:**

- **Sparse matrix functions**: The implementation uses `GridDef` objects instead of explicit bounds/resolution parameters for better type safety and consistency with the rest of the package.
- **Weight matrix generation**: The implementation accepts pre-processed weight columns and factors instead of dynamic filter functions, providing a more predictable API.
- **Output format**: The `writeEmis` function outputs structured DataFrames suitable for further processing, rather than directly writing shapefiles.

These differences make the implementation more composable and testable while maintaining full compatibility with the workflow described in the reference notebook. The core spatial processing workflow is complete and functional.
