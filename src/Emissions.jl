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

# Export all public functions and types
# Constants and unit conversions
export tonperyear, tonpermonth, foot, kelvin, Pollutants

# Data types
export EmissionsDataFrame, SurrogateSpec, GridDef, SpatialProcessor, Config, IndexInfo

# FF10 data formats
export FF10NonPointDataFrame, FF10PointDataFrame, FF10NonRoadDataFrame, FF10OnRoadDataFrame

# I/O functions
export strip_missing, getCountry, read_grid, getShapefilePath, validateShapefile, readSrgSpecSMOKE, NewSpatialProcessor

# Spatial processing
export NewPolygon, NewGridIrregular, setupSpatialProcessor, findCountyPolygon, GetIndex, recordToGrid, GridFactors, uniqueCoordinates, uniqueLoc

# Surrogate operations
export generate_data_sparse_matrices, generate_weight_sparse_matrices, generate_grid_sparse_matrices, generate_countySurrogate, update_locIndex

# Output functions
export writeEmis, find_surrogate_by_code, get_data_weight_shapefiles

end
