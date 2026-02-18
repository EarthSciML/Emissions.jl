# IOAPI/NetCDF Output

## Overview

Air quality models such as CMAQ and CAMx require input emissions in the IOAPI-3
(Input/Output Applications Programming Interface) NetCDF format. Emissions.jl can
write gridded emissions directly to IOAPI-format NetCDF files, enabling seamless
integration with these modeling systems.

```@docs
write_ioapi
```

## IOAPI Format Details

The output files follow the IOAPI-3 conventions:

- **Dimensions**: `TSTEP` (unlimited), `LAY`, `ROW`, `COL`, `VAR`, `DATE-TIME`
- **TFLAG variable**: Julian date (YYYYDDD) and time (HHMMSS) for each timestep
- **Pollutant variables**: One `Float32` variable per species with dimensions `(COL, ROW, LAY, TSTEP)`
- **Global attributes**: Grid projection parameters, domain extent, temporal metadata

## Usage

### Writing from Model-Ready Data

```@example ioapi
using Emissions, Dates

grid = NewGridRegular("TEST_GRID", 4, 3, "+proj=longlat +datum=WGS84",
    1.0, 1.0, -100.0, 30.0)

hours = [DateTime(2019, 7, 1, h) for h in 0:2]

model_data = Dict{String, Array{Float64, 4}}()
model_data["NO2"] = zeros(3, 4, 1, 3)
model_data["NO2"][1, 1, 1, 1] = 1.5
model_data["SO2"] = ones(3, 4, 1, 3) * 0.1

outfile = tempname() * ".nc"
write_ioapi(outfile, model_data, grid, hours;
    description="Example emissions", history="doc example")
println("Wrote IOAPI file: $outfile")
rm(outfile)
nothing # hide
```

### Writing from a Merged DataFrame

A convenience method accepts a merged `DataFrame` directly and converts it
using `to_model_ready` internally:

```@example ioapi
using DataFrames

merged = DataFrame(
    grid_row = [1, 1, 2],
    grid_col = [1, 2, 1],
    hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)],
    pollutant = ["CO", "CO", "CO"],
    emission_rate = [1.0, 2.0, 3.0]
)

outfile = tempname() * ".nc"
write_ioapi(outfile, merged, grid, [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)])
println("Wrote IOAPI file from DataFrame: $outfile")
rm(outfile)
nothing # hide
```

### Verifying Output

You can verify the output file using NCDatasets.jl:

```@example ioapi
using NCDatasets

# Write a test file
outfile = tempname() * ".nc"
model_data2 = Dict("CO" => ones(3, 4, 1, 1))
write_ioapi(outfile, model_data2, grid, [DateTime(2019, 7, 1, 0)])

ds = NCDatasets.Dataset(outfile, "r")
println("Grid type (GDTYP): ", ds.attrib["GDTYP"])
println("Grid name (GDNAM): ", strip(ds.attrib["GDNAM"]))
println("Dimensions: NCOLS=$(ds.attrib["NCOLS"]), NROWS=$(ds.attrib["NROWS"]), NLAYS=$(ds.attrib["NLAYS"])")
println("Variables: ", strip(ds.attrib["VAR-LIST"]))
close(ds)
rm(outfile)
nothing # hide
```
