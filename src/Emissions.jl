module Emissions
export g,
    ErrAboveModelTop, findLayer, calcDeltaH, ASME, calcDeltaHPrecomputed, ASMEPrecomputed

using CSV
using DataFrames
using Shapefile
using LibGEOS
using GeoInterface
using Proj
using SparseArrays
using Printf
using Unitful

include("plumerise.jl")
include("constants.jl")
include("types.jl")
include("ff10.jl")
include("io.jl")
include("spatial.jl")
include("surrogates.jl")
include("output.jl")
include("pipeline.jl")
include("workflow.jl")

# Export all public functions and types
# Constants and unit conversions
export tonperyear, tonpermonth, foot, kelvin, Pollutants

# Data types
export EmissionsDataFrame, SurrogateSpec, GridDef, SpatialProcessor, Config, IndexInfo

# FF10 data formats
export FF10NonPointDataFrame, FF10PointDataFrame, FF10NonRoadDataFrame, FF10OnRoadDataFrame

# I/O functions
export strip_missing, getCountry, normalize_country, read_grid, read_gridref,
    getShapefilePath, validateShapefile, readSrgSpecSMOKE, NewSpatialProcessor

# Spatial processing
export NewPolygon, NewGridIrregular, setupSpatialProcessor, findCountyPolygon, GetIndex, recordToGrid, GridFactors, uniqueCoordinates, uniqueLoc

# Surrogate operations
export generate_data_sparse_matrices, generate_weight_sparse_matrices, generate_grid_sparse_matrices, generate_countySurrogate, update_locIndex

# Output functions
export writeEmis, find_surrogate_by_code, get_data_weight_shapefiles

# Pipeline functions
export read_ff10, aggregate_emissions, filter_known_pollutants,
    map_pollutant_names!, assign_surrogates, build_data_weight_map

# Workflow orchestration
export location_key, compute_grid_indices, refine_indices_with_surrogates,
    allocate_emissions_to_grid, process_emissions_spatial

end
