# Spatial Allocation with ConservativeRegridding

## Overview

Emissions.jl provides spatial allocation using ConservativeRegridding.jl for
efficient, mathematically rigorous conservative regridding. This approach uses
precomputed intersection area matrices with spatial indexing (STRtree) for
better performance with large grids compared to manual polygon intersection loops.

```@docs
build_regridder
grid_polygons
generate_data_sparse_matrices_cr
generate_weight_sparse_matrices_cr
```

## Usage

### Building a Regridder

```julia
using Emissions
import GeoInterface as GI

# Define grid
grid = NewGridRegular("test", 10, 10, "EPSG:4326", 0.5, 0.5, -80.0, 35.0)

# Source geometries (e.g., county polygons)
source_geoms = [...]  # Vector of GeoInterface polygons

# Build reusable regridder
regridder = build_regridder(source_geoms, grid)
```

### Surrogate Generation

The `_cr` variants of the surrogate generation functions use ConservativeRegridding
internally for computing intersection areas:

```julia
# Using ConservativeRegridding (recommended for large grids)
data_matrices = generate_data_sparse_matrices_cr(
    "counties.shp", "FIPS", grid, "EPSG:4326"
)
weight_matrix = generate_weight_sparse_matrices_cr(
    "population.shp", ["POP"], [1.0], grid, "EPSG:4326"
)

# Original GeometryOps-based approach (compatible fallback)
data_matrices_orig = generate_data_sparse_matrices(
    "counties.shp", "FIPS", grid, "EPSG:4326"
)
```

## Performance

ConservativeRegridding.jl uses Sort-Tile-Recursive (STR) tree spatial indexing
to efficiently identify intersecting polygon pairs, avoiding the O(n*m) brute
force approach of checking every source geometry against every grid cell. This
provides significant speedups for large grids and many source polygons.
