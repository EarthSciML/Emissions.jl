# Inventory Validation

## Overview

The inventory validation module provides quality assurance checks for emissions
inventory data before processing. It detects common data issues such as duplicate
records, out-of-range values, and missing required fields.

This implements functionality analogous to SMOKE's inventory QA checks, ensuring
data integrity before the computationally expensive processing steps.

```@docs
ValidationResult
check_duplicates
check_ranges
check_completeness
validate_inventory
```

## Usage

### Basic Validation

```@example validation
using Emissions, DataFrames

# Create sample inventory data
emissions = DataFrame(
    FIPS = ["36001", "36001", "36005", "06001"],
    SCC = ["2103007000", "2103007000", "2103007001", "2103007000"],
    POLID = ["NOX", "NOX", "VOC", "SO2"],
    ANN_VALUE = [100.0, 50.0, 200.0, -10.0]
)

result = validate_inventory(emissions)
println("Valid: ", result.valid)
println("Warnings: ", result.warnings)
println("Errors: ", result.errors)
println("Duplicates found: ", result.n_duplicates)
println("Range issues: ", result.n_range_issues)
println("Missing fields: ", result.n_missing_fields)
```

### Individual Checks

```@example validation
# Check for duplicates
dups = check_duplicates(emissions)
println("Duplicate groups: ", nrow(dups))

# Check value ranges
range_issues = check_ranges(emissions)
println("Range issues: ", nrow(range_issues))

# Check completeness
missing_fields = check_completeness(emissions)
println("Records with missing fields: ", nrow(missing_fields))
```

### Pipeline Integration

Validation can be enabled in the processing pipeline via the `do_validate` keyword:

```julia
result = process_emissions(;
    inventory_files = [("inventory.csv", :nonpoint)],
    grid = my_grid,
    do_validate = true,  # Enable validation after aggregation
    # ... other options
)
```
