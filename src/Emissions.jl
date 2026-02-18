module Emissions
export g,
    ErrAboveModelTop, findLayer, calcDeltaH, ASME, calcDeltaHPrecomputed, ASMEPrecomputed

using CSV
import ConservativeRegridding
using DataFrames
using Dates
import GeoDataFrames
import GeoInterface as GI
import GeometryOps as GO
import NCDatasets
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
include("speciation.jl")
include("controls.jl")
include("elevpoint.jl")
include("laypoint.jl")
include("biogenic.jl")
include("validation.jl")
include("orl.jl")
include("pipeline.jl")
include("workflow.jl")
include("temporal.jl")
include("merge.jl")
include("ioapi.jl")
include("report.jl")

# Export all public functions and types
# Constants and unit conversions
export tonperyear, tonpermonth, foot, kelvin, Pollutants

# Data types
export EmissionsDataFrame, SurrogateSpec, GridDef, is_regular, SpatialProcessor, Config, IndexInfo

# FF10 data formats
export FF10NonPointDataFrame, FF10PointDataFrame, FF10NonRoadDataFrame, FF10OnRoadDataFrame

# ORL data formats
export ORLNonPointDataFrame, ORLPointDataFrame, ORLNonRoadDataFrame,
    ORLOnRoadDataFrame, ORLFireDataFrame

# I/O functions
export strip_missing, getCountry, normalize_country, read_grid, read_gridref,
    getShapefilePath, validateShapefile, readSrgSpecSMOKE, NewSpatialProcessor,
    read_griddesc, write_ioapi

# Spatial processing
export NewPolygon, NewGridIrregular, setupSpatialProcessor, findCountyPolygon, GetIndex,
    recordToGrid, GridFactors, uniqueCoordinates, uniqueLoc,
    cell_bounds, cell_polygon, cell_area,
    build_regridder, grid_polygons,
    new_polygon, new_grid_irregular, get_index, grid_factors

# Surrogate operations
export generate_data_sparse_matrices, generate_weight_sparse_matrices,
    generate_grid_sparse_matrices, generate_countySurrogate, update_locIndex,
    generate_data_sparse_matrices_cr, generate_weight_sparse_matrices_cr

# Output functions
export writeEmis, find_surrogate_by_code, get_data_weight_shapefiles

# Validation
export ValidationResult, check_duplicates, check_ranges, check_completeness, validate_inventory

# Pipeline functions
export read_ff10, read_orl, aggregate_emissions, filter_known_pollutants,
    map_pollutant_names!, assign_surrogates, build_data_weight_map,
    process_emissions

# Workflow orchestration
export location_key, compute_grid_indices, refine_indices_with_surrogates,
    allocate_emissions_to_grid, process_emissions_spatial

# Temporal allocation
export read_temporal_profiles, read_temporal_xref, temporal_allocate,
    read_day_specific, read_hour_specific

# Merge
export merge_emissions, merge_categories, merge_categories_tracked,
    merge_2d_3d, to_model_ready

# QA Reporting
export ReportConfig, summary_by_pollutant, summary_by_region, summary_by_scc,
    summary_by_time, compare_inventories, emissions_report

# Speciation
export read_gspro, read_gsref, build_speciation_matrix, speciate_emissions

# Controls
export ControlSpec, read_growth_factors, read_control_factors, apply_controls

# Elevated source identification
export ElevationCriteria, DEFAULT_ELEVATION_CRITERIA,
    analytical_plume_rise, classify_point_sources, group_stacks

# Vertical layer allocation
export LayerConfig, MetProfile, compute_layer_fractions,
    allocate_point_to_layers, laypoint

# Biogenic emissions
export BiogenicConfig, read_beld, read_emission_factors,
    temperature_adjustment, light_adjustment, compute_biogenic_emissions

end
