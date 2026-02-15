# Complete Tutorial

This tutorial demonstrates the complete workflow for processing EPA National Emissions Inventory (NEI) data using Emissions.jl, from raw FF10 files to gridded, spatially-allocated outputs.

## Overview

The emissions processing workflow consists of several key steps:
1. **Setup**: Configuration and data paths
2. **Load**: Read FF10 format emissions files
3. **Aggregate**: Combine and sum emissions by location and source
4. **Filter**: Keep only known pollutants and apply quality filters
5. **Spatial Setup**: Initialize spatial processing components
6. **Grid Reference**: Assign spatial surrogates to emissions
7. **Spatial Allocation**: Generate allocation matrices from shapefiles
8. **Output**: Write gridded emissions to shapefile

This tutorial shows two approaches: a **synthetic data example** for learning and testing, and instructions for processing the **full NEI dataset**.

## Synthetic Data Walkthrough

The following example demonstrates the complete workflow using synthetic data that can be run without downloading large datasets.

### Step 1: Create Synthetic Emissions Data

The first step is to create (or load) emissions data in the FF10 format used by the EPA.
FF10 nonpoint files have 45 columns including location identifiers (`FIPS` county codes),
source classification codes (`SCC`), pollutant IDs, and annual emission values in tons/year.

```@example complete_workflow
using Emissions
using DataFrames
using CSV
using Printf
using Unitful

synthetic_data = DataFrame(
    COUNTRY = ["0", "0", "0", "0"],
    FIPS = ["36001", "36001", "36005", "36005"],
    TRIBAL_CODE = ["0", "0", "0", "0"],
    CENSUS_TRACT = ["0", "0", "0", "0"],
    SHAPE_ID = ["0", "0", "0", "0"],
    SCC = ["2103007000", "2103007000", "2103007000", "2103007000"],
    EMIS_TYPE = ["", "", "", ""],
    POLID = ["NOX", "VOC", "NOX", "VOC"],
    ANN_VALUE = [150.5, 75.2, 200.1, 125.8],  # tons/year
    ANN_PCT_RED = [0.0, 0.0, 0.0, 0.0],
    CONTROL_IDS = ["", "", "", ""],
    CONTROL_MEASURES = ["", "", "", ""],
    CURRENT_COST = [0.0, 0.0, 0.0, 0.0],
    CUMULATIVE_COST = [0.0, 0.0, 0.0, 0.0],
    PROJECTION_FACTOR = [1.0, 1.0, 1.0, 1.0],
    REG_CODES = ["", "", "", ""],
    CALC_METHOD = [1, 1, 1, 1],
    CALC_YEAR = [2019, 2019, 2019, 2019],
    DATE_UPDATED = ["", "", "", ""],
    DATA_SET_ID = ["", "", "", ""],
    JAN_VALUE = [0.0, 0.0, 0.0, 0.0],
    FEB_VALUE = [0.0, 0.0, 0.0, 0.0],
    MAR_VALUE = [0.0, 0.0, 0.0, 0.0],
    APR_VALUE = [0.0, 0.0, 0.0, 0.0],
    MAY_VALUE = [0.0, 0.0, 0.0, 0.0],
    JUN_VALUE = [0.0, 0.0, 0.0, 0.0],
    JUL_VALUE = [0.0, 0.0, 0.0, 0.0],
    AUG_VALUE = [0.0, 0.0, 0.0, 0.0],
    SEP_VALUE = [0.0, 0.0, 0.0, 0.0],
    OCT_VALUE = [0.0, 0.0, 0.0, 0.0],
    NOV_VALUE = [0.0, 0.0, 0.0, 0.0],
    DEC_VALUE = [0.0, 0.0, 0.0, 0.0],
    JAN_PCTRED = [0.0, 0.0, 0.0, 0.0],
    FEB_PCTRED = [0.0, 0.0, 0.0, 0.0],
    MAR_PCTRED = [0.0, 0.0, 0.0, 0.0],
    APR_PCTRED = [0.0, 0.0, 0.0, 0.0],
    MAY_PCTRED = [0.0, 0.0, 0.0, 0.0],
    JUN_PCTRED = [0.0, 0.0, 0.0, 0.0],
    JUL_PCTRED = [0.0, 0.0, 0.0, 0.0],
    AUG_PCTRED = [0.0, 0.0, 0.0, 0.0],
    SEP_PCTRED = [0.0, 0.0, 0.0, 0.0],
    OCT_PCTRED = [0.0, 0.0, 0.0, 0.0],
    NOV_PCTRED = [0.0, 0.0, 0.0, 0.0],
    DEC_PCTRED = [0.0, 0.0, 0.0, 0.0],
    COMMENT = ["Synthetic example", "Synthetic example", "Synthetic example", "Synthetic example"]
)

first(synthetic_data[!, [:FIPS, :SCC, :POLID, :ANN_VALUE]], 4)
```

We have 4 records: NOX and VOC emissions from two New York counties (Albany=36001, Bronx=36005), all from the same source category (SCC 2103007000, residential heating).

### Step 2: Load Emissions with FF10 Format Validation

The `FF10NonPointDataFrame` constructor validates the data structure and automatically
converts emission values from tons/year to kg/s (SI units):

```@example complete_workflow
ff10_data = FF10NonPointDataFrame(synthetic_data)
processed_emis = ff10_data.df

DataFrame(
    Original_tons_yr = synthetic_data.ANN_VALUE,
    Converted_kg_s = processed_emis.ANN_VALUE,
    Ratio = processed_emis.ANN_VALUE ./ synthetic_data.ANN_VALUE
)
```

The conversion factor is consistent across all records, as expected for a simple unit conversion from tons/year to kg/s.

### Step 3: Aggregate and Filter Emissions

Next we group emissions by key fields and sum duplicate records, then filter to keep only
known pollutants and map them to standardized names:

```@example complete_workflow
grouped_emis = combine(
    groupby(processed_emis, [:POLID, :COUNTRY, :FIPS, :SCC]),
    :ANN_VALUE => sum => :ANN_VALUE
)

known_polls = filter(row -> haskey(Pollutants, row.POLID), grouped_emis)
known_polls[!, :POLID] = [Pollutants[p] for p in known_polls[!, :POLID]]

known_polls
```

The `Pollutants` dictionary maps FF10 pollutant codes (e.g., `"NOX"`, `"VOC"`) to
standardized names used downstream.

### Step 4: Set Up Spatial Processing Configuration

Before we can assign emissions to grid cells, we need supporting data structures: a
grid reference table that maps each `(COUNTRY, FIPS, SCC)` combination to a spatial
surrogate code, a surrogate specification describing how to spatially distribute emissions,
and a configuration object with file paths:

```@example complete_workflow
grid_ref = DataFrame(
    COUNTRY = ["USA", "USA", "USA"],
    FIPS = ["36001", "36005", "00000"],  # "00000" is a fallback surrogate
    SCC = ["2103007000", "2103007000", "2103007000"],
    Surrogate = [100, 100, 100]
)

grid_ref
```

```@example complete_workflow
srg_spec = SurrogateSpec(
    "USA",                    # Region
    "Population",             # Name
    100,                      # Code
    "population.shp",         # DataShapefile
    "POP2019",               # DataAttribute
    "area.shp",              # WeightShapefile
    "Area-based allocation",  # Details
    String[],                # BackupSurrogateNames
    ["AREA"],                # WeightColumns
    [1.0],                   # WeightFactors
    "",                      # FilterFunction
    String[],                # MergeNames
    Float64[]                # MergeMultipliers
)

println("Surrogate $(srg_spec.Code): $(srg_spec.Name)")
println("  Data shapefile: $(srg_spec.DataShapefile) [$(srg_spec.DataAttribute)]")
println("  Weight shapefile: $(srg_spec.WeightShapefile) [$(join(srg_spec.WeightColumns, ", "))]")
```

```@example complete_workflow
config = Config(
    ["/tmp/synthetic_gridref.csv"],      # f_gridRef
    "/tmp/synthetic_srgspec.csv",        # SrgSpec
    "/tmp/shapefiles/",                  # SrgShapefileDirectory
    "+proj=longlat +datum=WGS84",        # InputSR
    "+proj=lcc +lat_1=33 +lat_2=45",     # OutputSR
    "/tmp/synthetic_grid.txt",           # GridFile
    "SyntheticGrid",                     # GridName
    "/tmp/counties.shp",                 # Counties
    "/tmp/output/"                       # EmisShp
)

println("Grid name: $(config.GridName)")
println("Input projection: $(config.InputSR)")
println("Output projection: $(config.OutputSR)")
```

### Step 5: Assign Spatial Surrogates

We join the emissions data with the grid reference to assign each record a spatial
surrogate code. The surrogate determines how emissions from a county are distributed
across grid cells (e.g., proportional to population):

```@example complete_workflow
emissions_with_surrogates = leftjoin(
    known_polls,
    grid_ref,
    on = [:COUNTRY, :FIPS, :SCC]
)

emissions_with_surrogates
```

Any rows with `missing` in the `Surrogate` column would need fallback handling. In real
processing, these would be matched against `FIPS="00000"` records which provide state-level
default surrogates.

```@example complete_workflow
unmatched = filter(row -> ismissing(row.Surrogate), emissions_with_surrogates)
println("Unmatched emissions requiring fallback: $(nrow(unmatched))")
```

### Step 6: Key Constants and Unit Conversions

Emissions.jl provides built-in unit conversion constants and a temperature conversion
function used throughout the processing pipeline:

```@example complete_workflow
DataFrame(
    Conversion = ["tons/year to kg/s", "tons/month to kg/s", "feet to meters"],
    Factor = [ustrip(tonperyear), ustrip(tonpermonth), ustrip(foot)]
)
```

```@example complete_workflow
temps_f = [32.0, 68.0, 212.0]
DataFrame(
    Fahrenheit = temps_f,
    Kelvin = [ustrip(kelvin(t)) for t in temps_f],
    Description = ["Freezing point", "Room temperature", "Boiling point"]
)
```

The `Pollutants` dictionary maps all supported FF10 pollutant codes to their standardized names:

```@example complete_workflow
DataFrame(
    FF10_Name = collect(keys(Pollutants)),
    Standard_Name = collect(values(Pollutants))
)
```

### Step 7: Spatial Processing Overview

After surrogate assignment, the remaining spatial processing steps distribute emissions
to grid cells using spatial surrogate data derived from shapefiles. Here is a summary of
the workflow and its key functions:

| Step | Function | Description |
|:-----|:---------|:------------|
| Load shapefiles | `generate_data_sparse_matrices()` | Read population/area data |
| Weight shapefiles | `generate_weight_sparse_matrices()` | Read surrogate weights |
| Grid matrices | `generate_grid_sparse_matrices()` | Create target grid coverage |
| Surrogate generation | `generate_countySurrogate()` | Compute normalized allocation fractions |
| Location indexing | `GetIndex()` | Find grid cells for each emission location |
| Final allocation | `recordToGrid()` | Distribute emissions to cells |
| Output | `writeEmis()` | Write final gridded emissions |

An example of what the final gridded output looks like:

```@example complete_workflow
DataFrame(
    GridCell = [1, 1, 2, 2],
    Pollutant = ["NOX", "VOC", "NOX", "VOC"],
    EmissionRate_kg_s = [1.2e-6, 6.1e-7, 1.6e-6, 1.0e-6],
    SourceCount = [2, 2, 2, 2],
    County = ["36001", "36001", "36005", "36005"]
)
```

For a fully executable example of spatial allocation with synthetic surrogates, see the
[Spatial Processing Pipeline](@ref) page.

## Processing the Full NEI Dataset

For processing the complete EPA National Emissions Inventory, follow these steps:

### Step 1: Download NEI Data

```julia
# The full NEI dataset can be downloaded from EPA's website:
# https://www.epa.gov/air-emissions-inventories/2019-national-emissions-inventory-nei-data

# Example file structure for 2019 NEI:
nei_base_url = "https://gaftp.epa.gov/air/nei/2019/"
file_downloads = [
    "nonpoint/2019v1_nonpoint_20210121_csv.zip",
    "point/2019v1_point_20210121_csv.zip",
    "nonroad/2019v1_nonroad_20210121_csv.zip",
    "onroad/2019v1_onroad_20210121_csv.zip"
]

# Download and extract files
using HTTP, ZipFile
for file in file_downloads
    url = nei_base_url * file
    local_path = joinpath("data", "nei2019", basename(file))

    # Download
    println("Downloading: $url")
    HTTP.download(url, local_path)

    # Extract ZIP file
    zip_reader = ZipFile.Reader(local_path)
    for file_info in zip_reader.files
        # Extract CSV files to data directory
        if endswith(file_info.name, ".csv")
            extracted_path = joinpath("data", "nei2019", file_info.name)
            open(extracted_path, "w") do io
                write(io, read(file_info))
            end
            println("  Extracted: $(file_info.name)")
        end
    end
    close(zip_reader)
end
```

### Step 2: Set Up Real Configuration

```julia
# Configure paths for full NEI processing
real_config = Config(
    # Grid reference files (multiple regions)
    [
        "data/spatial/gridref_usa_2019.csv",
        "data/spatial/gridref_canada_2019.csv",
        "data/spatial/gridref_mexico_2019.csv"
    ],

    # Surrogate specification file
    "data/spatial/surrogate_specification_2019.csv",

    # Directory containing all spatial surrogate shapefiles
    "data/spatial/shapefiles/",

    # Input coordinate system (geographic coordinates)
    "+proj=longlat +datum=WGS84 +no_defs",

    # Output coordinate system (Lambert Conformal Conic for US modeling)
    "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",

    # Grid definition file (defines target modeling grid)
    "data/grids/inmap_conus_grid.txt",

    # Grid identifier
    "InMAP_CONUS",

    # County boundaries shapefile
    "data/spatial/shapefiles/tl_2019_us_county.shp",

    # Output directory for gridded emissions
    "output/nei2019_gridded/"
)
```

### Step 3: Process Full Dataset by Sector

```julia
# Define readers for each emission sector
readers = Dict(
    "nonpoint" => (f -> FF10NonPointDataFrame(CSV.read(f, DataFrame, comment="#"))),
    "point" => (f -> FF10PointDataFrame(CSV.read(f, DataFrame, comment="#"))),
    "nonroad" => (f -> FF10NonRoadDataFrame(CSV.read(f, DataFrame, comment="#"))),
    "onroad" => (f -> FF10OnRoadDataFrame(CSV.read(f, DataFrame, comment="#")))
)

# Process each sector separately to manage memory
all_emissions = DataFrame[]

for (sector, reader) in readers
    println("Processing $sector emissions...")

    # Find all CSV files for this sector
    sector_files = filter(f -> contains(f, sector) && endswith(f, ".csv"),
                         readdir("data/nei2019/", join=true))

    sector_emissions = DataFrame[]
    for file in sector_files
        println("  Reading: $(basename(file))")
        try
            ff10_data = reader(file)
            push!(sector_emissions, ff10_data.df)
        catch e
            println("    Warning: Error reading $file: $e")
        end
    end

    # Concatenate all files in this sector
    if !isempty(sector_emissions)
        sector_combined = vcat(sector_emissions...)
        println("  $sector total: $(nrow(sector_combined)) records")
        push!(all_emissions, sector_combined)
    end
end

# Combine all sectors
println("Combining all emission sectors...")
combined_emissions = vcat(all_emissions...)
println("Total emissions: $(nrow(combined_emissions)) records")
```

### Step 4: Full Spatial Processing

```julia
# Initialize spatial processing for the full domain
sp, gd = setupSpatialProcessor(real_config)

# Aggregate emissions by location and source characteristics
aggregated = combine(
    groupby(combined_emissions, [:POLID, :COUNTRY, :FIPS, :SCC, :LON, :LAT]),
    :ANN_VALUE => sum => :ANN_VALUE,
    :STACK_HEIGHT => mean => :STACK_HEIGHT,  # For point sources
    :STACK_DIAMETER => mean => :STACK_DIAMETER,
    :STACK_TEMP => mean => :STACK_TEMP,
    :STACK_FLOW => sum => :STACK_FLOW,
    :STACK_VELOCITY => mean => :STACK_VELOCITY
)

# Filter to known pollutants and standardize names
known_emis = filter(row -> haskey(Pollutants, row.POLID), aggregated)
known_emis[!, :POLID] = [Pollutants[p] for p in known_emis[!, :POLID]]

# Join with grid reference to assign surrogates
emissions_with_surrogates = leftjoin(known_emis, sp.GridRef, on=[:COUNTRY, :FIPS, :SCC])

# Handle fallback surrogates for unmatched emissions
unmatched = filter(row -> ismissing(row.Surrogate), emissions_with_surrogates)
if nrow(unmatched) > 0
    # Try fallback with FIPS = "00000" (statewide surrogates)
    fallback_ref = filter(row -> row.FIPS == "00000", sp.GridRef)

    # Match unmatched emissions with fallback surrogates
    unmatched_fallback = leftjoin(
        select(unmatched, Not(:Surrogate)),
        select(fallback_ref, [:COUNTRY, :SCC, :Surrogate => :FallbackSurrogate]),
        on = [:COUNTRY, :SCC]
    )
    unmatched_fallback[!, :Surrogate] = unmatched_fallback[!, :FallbackSurrogate]
    select!(unmatched_fallback, Not(:FallbackSurrogate))

    # Combine matched and fallback-matched emissions
    matched = filter(row -> !ismissing(row.Surrogate), emissions_with_surrogates)
    final_emissions = vcat(matched, unmatched_fallback)
else
    final_emissions = emissions_with_surrogates
end

println("Final emissions ready for spatial allocation: $(nrow(final_emissions)) records")
```

### Step 5: Generate Full Spatial Allocation

```julia
# Define processing domain bounds (full CONUS)
conus_bounds = (
    x_min = -130.0,  # Western boundary
    x_max = -65.0,   # Eastern boundary
    y_min = 20.0,    # Southern boundary
    y_max = 55.0     # Northern boundary
)

resolution = 0.01  # ~1km resolution

# Generate sparse matrices for all shapefiles used in surrogates
println("Generating spatial allocation matrices...")

all_sparse_data = Dict()
all_sparse_weight = Dict()

for spec in sp.SrgSpecs
    println("Processing surrogate $(spec.Code): $(spec.Name)")

    # Generate data matrices
    if !haskey(all_sparse_data, spec.DataShapefile)
        data_matrices = generate_data_sparse_matrices(
            conus_bounds.x_min, conus_bounds.x_max,
            conus_bounds.y_min, conus_bounds.y_max,
            resolution,
            joinpath(real_config.SrgShapefileDirectory, spec.DataShapefile),
            spec.DataAttribute
        )
        all_sparse_data[spec.DataShapefile] = data_matrices
    end

    # Generate weight matrices
    if !haskey(all_sparse_weight, spec.WeightShapefile)
        weight_matrices = generate_weight_sparse_matrices(
            conus_bounds.x_min, conus_bounds.x_max,
            conus_bounds.y_min, conus_bounds.y_max,
            resolution,
            joinpath(real_config.SrgShapefileDirectory, spec.WeightShapefile),
            spec.WeightColumns[1];  # Use first weight column
            Filter = spec.FilterFunction
        )
        all_sparse_weight[spec.WeightShapefile] = weight_matrices
    end
end
```

### Step 6: Write Final Output

```julia
# Generate county surrogates for all combinations
county_surrogates = Dict()
for spec in sp.SrgSpecs
    key = (spec.DataShapefile, spec.WeightShapefile)
    if !haskey(county_surrogates, key)
        county_surrogates[key] = generate_countySurrogate(
            all_sparse_data[spec.DataShapefile],
            all_sparse_weight[spec.WeightShapefile]
        )
    end
end

# For each unique emission location, create spatial index
unique_locations = unique(select(final_emissions, [:COUNTRY, :FIPS, :LON, :LAT]))
location_indices = Dict()

for row in eachrow(unique_locations)
    if !ismissing(row.LON) && !ismissing(row.LAT)
        # Point source - use coordinates
        location_key = (row.COUNTRY, row.FIPS, row.LON, row.LAT)
        location_indices[location_key] = GetIndex(row.LON, row.LAT, gd)
    else
        # Area source - use county polygon
        location_key = (row.COUNTRY, row.FIPS, "COUNTY")
        county_poly = findCountyPolygon(real_config.Counties, row.FIPS)
        if county_poly !== nothing
            location_indices[location_key] = GetIndex(county_poly, gd)
        end
    end
end

# Apply spatial allocation and write output
println("Writing final gridded emissions...")

output_records = []
for row in eachrow(final_emissions)
    if !ismissing(row.Surrogate)
        # Find appropriate surrogate data
        spec = find_surrogate_by_code(sp.SrgSpecs, row.Surrogate)
        if spec !== nothing
            # Get location index
            if !ismissing(row.LON) && !ismissing(row.LAT)
                location_key = (row.COUNTRY, row.FIPS, row.LON, row.LAT)
            else
                location_key = (row.COUNTRY, row.FIPS, "COUNTY")
            end

            if haskey(location_indices, location_key)
                index_info = location_indices[location_key]

                # Apply surrogate allocation
                surrogate_key = (spec.DataShapefile, spec.WeightShapefile)
                if haskey(county_surrogates, surrogate_key)
                    surrogate_data = county_surrogates[surrogate_key]

                    # Distribute emissions to grid cells
                    allocated = recordToGrid(
                        row.ANN_VALUE,
                        index_info,
                        get(surrogate_data, row.FIPS, nothing)
                    )

                    # Store results for output
                    for (cell_idx, emission_amount) in pairs(allocated.nzval)
                        push!(output_records, (
                            grid_cell = cell_idx,
                            pollutant = row.POLID,
                            emission_kg_s = emission_amount,
                            source_fips = row.FIPS,
                            source_scc = row.SCC
                        ))
                    end
                end
            end
        end
    end
end

# Convert to DataFrame and write shapefile
output_df = DataFrame(output_records)
writeEmis(
    joinpath(real_config.EmisShp, "nei2019_gridded_emissions.shp"),
    gd,
    output_df,
    real_config.OutputSR,
    "Final gridded 2019 NEI emissions"
)

println("Processing complete!")
println("Output written to: $(real_config.EmisShp)")
println("Total gridded emission records: $(nrow(output_df))")
```

## Performance Considerations

For processing the full NEI dataset (~50GB total):

- **Memory**: 32-64GB RAM recommended for full domain processing
- **Storage**: ~200GB for input data, intermediate files, and outputs
- **Processing Time**: 4-8 hours for complete CONUS domain at 1km resolution
- **Parallelization**: Consider geographic tiling or surrogate-based parallelization

For smaller domains or testing, use geographic bounds to subset the processing area and reduce computational requirements.

## Key Functions Reference

### FF10 Data Loading
- `FF10NonPointDataFrame()` - Area source emissions (45 columns)
- `FF10PointDataFrame()` - Point source emissions (77 columns)
- `FF10NonRoadDataFrame()` - Non-road mobile emissions (45 columns)
- `FF10OnRoadDataFrame()` - On-road mobile emissions (45 columns)

### Configuration and Setup
- `Config()` - Main configuration structure
- `setupSpatialProcessor()` - Initialize spatial processing

### Spatial Allocation
- `generate_data_sparse_matrices()` - Load data shapefiles
- `generate_weight_sparse_matrices()` - Load weight shapefiles
- `generate_countySurrogate()` - Create allocation surrogates
- `GetIndex()` - Map locations to grid cells
- `recordToGrid()` - Distribute emissions to cells

### Output and Utilities
- `writeEmis()` - Write final gridded shapefile
- `find_surrogate_by_code()` - Look up surrogate specifications
- Unit conversion constants: `tonperyear`, `tonpermonth`, `foot`, `kelvin()`
- Pollutant mapping: `Pollutants` dictionary