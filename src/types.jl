export EmissionsDataFrame,
    SurrogateSpec, GridDef, SpatialProcessor, Config, IndexInfo

"""
    EmissionsDataFrame

Abstract type for emissions data frames wrapping FF10 format data.
"""
abstract type EmissionsDataFrame end

"""
    SurrogateSpec

Holds surrogate specification data for spatial allocation of emissions.
"""
struct SurrogateSpec
    Region::String
    Name::String
    Code::Int
    DataShapefile::String
    DataAttribute::String
    WeightShapefile::String
    Details::String
    BackupSurrogateNames::Vector{String}
    WeightColumns::Vector{String}
    WeightFactors::Vector{Float64}
    FilterFunction::String
    MergeNames::Vector{String}
    MergeMultipliers::Vector{Float64}
end

"""
    GridDef

Specifies the grid that we are allocating the emissions to.
Each cell is described by its bounding box in `Extent`:
`[(xmin, ymin), (xmax, ymax)]`.
"""
struct GridDef
    Name::String
    Nx::Int
    Ny::Int
    SR::String
    Extent::Vector{Vector{Tuple{Float64, Float64}}}
end

"""
    SpatialProcessor

Spatializes emissions records using surrogate specifications and grid definitions.
"""
struct SpatialProcessor
    SrgSpecs::Vector{SurrogateSpec}
    Grids::GridDef
    GridRef::DataFrame
    InputSR::String
    MatchFullSCC::Bool
    MemCacheSize::Int
    MaxMergeDepth::Int
end

"""
    Config

Holds configuration data for the emissions processing pipeline.

# Fields
- `f_gridRef`: Paths to grid reference files
- `SrgSpec`: Path to surrogate specification file
- `SrgShapefileDirectory`: Directory containing surrogate shapefiles
- `InputSR`: Input spatial reference (projection string)
- `OutputSR`: Output spatial reference (projection string)
- `GridFile`: Path to grid definition file
- `GridName`: Name of the grid
- `Counties`: Path to counties shapefile
- `EmisShp`: Path to output emissions shapefile directory
"""
struct Config
    f_gridRef::Vector{String}
    SrgSpec::String
    SrgShapefileDirectory::String
    InputSR::String
    OutputSR::String
    GridFile::String
    GridName::String
    Counties::String
    EmisShp::String
end

"""
    IndexInfo

Holds grid index information for gridded emissions, including which grid cells
a source maps to and the fraction of emissions allocated to each cell.
"""
struct IndexInfo
    rows::Vector{Int}
    cols::Vector{Int}
    fracs::Vector{Float64}
    inGrid::Bool
    coveredByGrid::Bool
end
