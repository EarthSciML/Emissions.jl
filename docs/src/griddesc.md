# GRIDDESC File Support

## Overview

The GRIDDESC (Grid Description) file is the standard IOAPI format for defining model grids
used by air quality models such as CMAQ, CAMx, and SMOKE. Emissions.jl can read GRIDDESC
files and create `GridDef` objects for spatial allocation of emissions.

```@docs
read_griddesc
```

## GRIDDESC File Format

A GRIDDESC file has two segments separated by delimiter lines (a line containing only `' '`):

**Segment 1** defines coordinate systems. Each entry consists of two lines:
1. Coordinate system name (quoted string)
2. Six numeric parameters: `coordtype P_ALP P_BET P_GAM XCENT YCENT`

**Segment 2** defines grids. Each entry consists of two lines:
1. Grid name (quoted string)
2. Eight parameters: `coord_name XORIG YORIG XCELL YCELL NCOLS NROWS NTHIK`

### Supported Coordinate Types

| Code | Name     | Description                    |
|------|----------|--------------------------------|
| 1    | LATGRD3  | Geographic latitude-longitude  |
| 2    | LAMGRD3  | Lambert Conformal Conic        |
| 5    | UTMGRD3  | Universal Transverse Mercator  |
| 6    | POLGRD3  | Polar Stereographic            |
| 7    | EQMGRD3  | Equatorial Mercator            |
| 8    | TRMGRD3  | Transverse Mercator            |

## Usage

### Reading a GRIDDESC File

```@example griddesc
using Emissions

# Create an example GRIDDESC file
griddesc_content = """
' '
'LCC_US'
2 33.0 45.0 -97.0 -97.0 40.0
' '
'US36_CRO'
'LCC_US' -2736000.0 -2088000.0 36000.0 36000.0 148 112 1
' '
"""

path = tempname()
write(path, griddesc_content)

grid = read_griddesc(path, "US36_CRO")
println("Grid: $(grid.Name)")
println("Size: $(grid.Nx) x $(grid.Ny)")
println("Cell size: $(grid.Dx) x $(grid.Dy)")
println("Origin: ($(grid.X0), $(grid.Y0))")
println("Projection: $(grid.SR)")
rm(path)
nothing # hide
```

### Using with Spatial Processing

The `setupSpatialProcessor` function automatically detects and reads GRIDDESC files.
When a grid file is provided in the `Config`, it first attempts to parse it as a GRIDDESC
file before falling back to other formats.
