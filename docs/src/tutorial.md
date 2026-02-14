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

```@example complete_workflow
using Emissions
using DataFrames
using CSV
using Printf
using Unitful

# =============================================================================
# Step 1: Create Synthetic Emissions Data
# =============================================================================
println("=== Step 1: Creating Synthetic Emissions Data ===")

# Create synthetic data matching exact FF10 nonpoint format (45 columns)
synthetic_data = DataFrame(
    COUNTRY = ["0", "0", "0", "0"],
    FIPS = ["36001", "36001", "36005", "36005"],
    TRIBAL_CODE = ["0", "0", "0", "0"],
    CENSUS_TRACT = ["0", "0", "0", "0"],
    SHAPE_ID = ["0", "0", "0", "0"],
    SCC = ["2103007000", "2103007000", "2103007000", "2103007000"],
    EMIS_TYPE = ["", "", "", ""],
    POLID = ["NOX", "VOC", "NOX", "VOC"],                       # Pollutant ID
    ANN_VALUE = [150.5, 75.2, 200.1, 125.8],                   # tons/year
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

println("Created synthetic emissions with $(nrow(synthetic_data)) records and $(ncol(synthetic_data)) columns")
println("First few rows:")
show(first(synthetic_data[!, 1:8], 4), allcols=true)  # Show first 8 columns

# =============================================================================
# Step 2: Load Emissions with FF10 Format Validation
# =============================================================================
println("\n\n=== Step 2: Loading Emissions with FF10 Format Validation ===")

# The FF10NonPointDataFrame constructor validates structure and converts units
ff10_data = FF10NonPointDataFrame(synthetic_data)
processed_emis = ff10_data.df
println("✅ Successfully loaded and validated FF10 nonpoint data")
println("Shape: $(nrow(processed_emis)) rows × $(ncol(processed_emis)) columns")

# Show the unit conversion in action
println("\nUnit conversion examples:")
println("Original tons/year: $(synthetic_data[1, :ANN_VALUE])")
println("Converted to kg/s: $(processed_emis[1, :ANN_VALUE])")
@printf "Conversion factor: %.2e kg/s per ton/year\n" (processed_emis[1, :ANN_VALUE] / synthetic_data[1, :ANN_VALUE])

# =============================================================================
# Step 3: Aggregating and Filtering Emissions
# =============================================================================
println("\n\n=== Step 3: Aggregating and Filtering Emissions ===")

# Group by key fields and sum emissions (this removes duplicates)
grouped_emis = combine(
    groupby(processed_emis, [:POLID, :COUNTRY, :FIPS, :SCC]),
    :ANN_VALUE => sum => :ANN_VALUE
)

println("After grouping: $(nrow(grouped_emis)) unique emission records")
show(grouped_emis, allcols=true)

# Filter to keep only known pollutants
known_polls = filter(row -> haskey(Pollutants, row.POLID), grouped_emis)
println("\n\nKnown pollutants found: $(nrow(known_polls)) records")

# Map to standard pollutant names
known_polls[!, :POLID] = [Pollutants[p] for p in known_polls[!, :POLID]]
println("Mapped to standard names:")
show(known_polls, allcols=true)

# =============================================================================
# Step 4: Create Supporting Data Structures
# =============================================================================
println("\n\n=== Step 4: Setting Up Spatial Processing Configuration ===")

# Create synthetic grid reference data
grid_ref = DataFrame(
    COUNTRY = ["USA", "USA", "USA"],
    FIPS = ["36001", "36005", "00000"],  # Include fallback (00000) surrogate
    SCC = ["2103007000", "2103007000", "2103007000"],
    Surrogate = [100, 100, 100]  # All use surrogate code 100
)

# Create synthetic surrogate specification
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

# Create synthetic configuration (in real use, these would be actual file paths)
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

println("Grid Reference data:")
show(grid_ref, allcols=true)

println("\n\nSurrogate Specifications:")
println("Surrogate $(srg_spec.Code): $(srg_spec.Name)")
println("  Data: $(srg_spec.DataShapefile) [$(srg_spec.DataAttribute)]")
println("  Weight: $(srg_spec.WeightShapefile) [$(join(srg_spec.WeightColumns, ", "))]")
println("  Details: $(srg_spec.Details)")

println("\n\nConfiguration:")
println("Input projection: $(config.InputSR)")
println("Output projection: $(config.OutputSR)")
println("Grid name: $(config.GridName)")

# =============================================================================
# Step 5: Assign Spatial Surrogates
# =============================================================================
println("\n\n=== Step 5: Assigning Spatial Surrogates ===")

# POLID column already matches grid reference format (no renaming needed)

# Join with grid reference to assign surrogates
emissions_with_surrogates = leftjoin(
    known_polls,
    grid_ref,
    on = [:COUNTRY, :FIPS, :SCC]
)

println("Emissions with assigned surrogates:")
show(emissions_with_surrogates, allcols=true)

# Check for unmatched emissions (would need fallback surrogates in real data)
unmatched = filter(row -> ismissing(row.Surrogate), emissions_with_surrogates)
println("\nUnmatched emissions requiring fallback surrogates: $(nrow(unmatched))")

if nrow(unmatched) > 0
    println("In real processing, these would be matched against FIPS='00000' records")
end

# =============================================================================
# Step 6: Demonstrate Key Constants and Unit Conversions
# =============================================================================
println("\n\n=== Step 6: Key Constants and Unit Conversions ===")

println("Emissions unit conversion factors:")
@printf "tons/year → kg/s: %.2e\n" ustrip(tonperyear)
@printf "tons/month → kg/s: %.2e\n" ustrip(tonpermonth)
@printf "feet → meters: %.4f\n" ustrip(foot)

println("\nTemperature conversion examples:")
temps_f = [32.0, 68.0, 212.0]  # Freezing, room temp, boiling
for temp_f in temps_f
    temp_k = kelvin(temp_f)
    @printf "%.1f°F = %.2f K\n" temp_f ustrip(temp_k)
end

println("\nSupported pollutants:")
for (ff10_name, standard_name) in pairs(Pollutants)
    println("  $ff10_name → $standard_name")
end

# =============================================================================
# Step 7: Spatial Processing Workflow Overview
# =============================================================================
println("\n\n=== Step 7: Spatial Processing Workflow (Conceptual) ===")

println("""
The complete spatial processing workflow would continue with:

1. **Shapefile Loading**: Read population and area shapefiles
   - generate_data_sparse_matrices() for population data
   - generate_weight_sparse_matrices() for area weights
   - generate_grid_sparse_matrices() for target grid

2. **Surrogate Generation**:
   - generate_countySurrogate() combines data and weights
   - Creates normalized allocation fractions that sum to 1.0

3. **Location Indexing**:
   - GetIndex() finds grid cells for each emission location
   - Handles partial grid cell coverage with fractions

4. **Final Allocation**:
   - recordToGrid() distributes emissions to grid cells
   - Applies surrogate fractions and location fractions

5. **Output Generation**:
   - writeEmis() creates final gridded shapefile
   - Includes emission amounts, grid coordinates, metadata
""")

# Show what the final output structure would look like
final_output_structure = DataFrame(
    GridCell = [1, 1, 2, 2],
    Pollutant = ["NOX", "VOC", "NOX", "VOC"],
    EmissionRate_kg_s = [1.2e-6, 6.1e-7, 1.6e-6, 1.0e-6],
    SourceCount = [2, 2, 2, 2],
    County = ["36001", "36001", "36005", "36005"]
)

println("Final gridded output structure:")
show(final_output_structure, allcols=true)

println("\n✅ Synthetic workflow demonstration complete!")
```

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