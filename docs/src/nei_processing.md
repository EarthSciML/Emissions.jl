# NEI Processing

## Overview

The Emissions.jl package provides tools for processing EPA National Emissions Inventory (NEI)
data. It supports reading FF10-format emissions files, spatial allocation using surrogate
shapefiles, and gridding emissions to model-ready formats.

The package also includes plume rise calculations based on ASME (1973), as described in
Seinfeld and Pandis, "Atmospheric Chemistry and Physics - From Air Pollution to Climate Change".

## Plume Rise

The following functions implement the ASME (1973) plume rise algorithm:

- `findLayer` - Find the model layer containing a given height
- `calcDeltaH` - Calculate plume rise
- `ASME` - Calculate effective emissions height with plume rise
- `calcDeltaHPrecomputed` - Calculate plume rise with precomputed meteorological parameters
- `ASMEPrecomputed` - Calculate effective emissions height with precomputed parameters

## Constants and Unit Conversions

```@docs
tonperyear
tonpermonth
foot
kelvin
Pollutants
```

## Data Types

```@docs
EmissionsDataFrame
SurrogateSpec
GridDef
SpatialProcessor
Config
IndexInfo
```

## FF10 Data Formats

The EPA FF10 (Flat File 10) format is the standard format for emissions inventory data.
The package supports four FF10 format types:

```@docs
FF10NonPointDataFrame
FF10PointDataFrame
FF10NonRoadDataFrame
FF10OnRoadDataFrame
```

## I/O Functions

```@docs
strip_missing
getCountry
read_grid
getShapefilePath
validateShapefile
readSrgSpecSMOKE
NewSpatialProcessor
```

## Spatial Processing

```@docs
NewPolygon
NewGridRegular
NewGridIrregular
setupSpatialProcessor
findCountyPolygon
GetIndex
recordToGrid
GridFactors
uniqueCoordinates
uniqueLoc
```

## Surrogate Operations

```@docs
generate_data_sparse_matrices
generate_weight_sparse_matrices
generate_grid_sparse_matrices
generate_countySurrogate
update_locIndex
```

## Output

```@docs
find_surrogate_by_code
get_data_weight_shapefiles
writeEmis
```

## Pipeline Functions

See the [Spatial Processing Pipeline](@ref) page for documentation on the mid-level pipeline
functions: `read_ff10`, `normalize_country`, `read_gridref`, `aggregate_emissions`,
`filter_known_pollutants`, `map_pollutant_names!`, `assign_surrogates`, and `build_data_weight_map`.
