export EmissionsDataFrame,
    SurrogateSpec, GridDef, is_regular, SpatialProcessor, Config, IndexInfo

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

For regular grids (uniform dx, dy), optional fields `Dx`, `Dy`, `X0`, `Y0`
enable O(1) point-in-cell lookup instead of linear scan.
"""
struct GridDef
    Name::String
    Nx::Int
    Ny::Int
    SR::String
    Extent::Vector{Vector{Tuple{Float64, Float64}}}
    Dx::Float64
    Dy::Float64
    X0::Float64
    Y0::Float64

    # Regular grid constructor (with dx, dy, x0, y0)
    function GridDef(name, nx, ny, sr, extent, dx, dy, x0, y0)
        return new(name, nx, ny, sr, extent, dx, dy, x0, y0)
    end

    # Irregular grid constructor (no dx, dy info - uses NaN sentinel)
    function GridDef(name, nx, ny, sr, extent)
        return new(name, nx, ny, sr, extent, NaN, NaN, NaN, NaN)
    end
end

"""
    is_regular(grid::GridDef) -> Bool

Check if a grid has regular spacing (uniform dx, dy).
"""
is_regular(grid::GridDef) = !isnan(grid.Dx)

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
