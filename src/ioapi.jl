export write_ioapi

"""
    _parse_proj4_to_ioapi(sr::AbstractString) -> NamedTuple

Parse a Proj4 projection string and return IOAPI-compatible global attributes
as a NamedTuple with fields `GDTYP`, `P_ALP`, `P_BET`, `P_GAM`, `XCENT`, `YCENT`.

This is the inverse of `_coordtype_to_proj4`.
"""
function _parse_proj4_to_ioapi(sr::AbstractString)
    # Extract projection type
    m = match(r"\+proj=(\w+)", sr)
    proj = m === nothing ? "" : m.captures[1]

    # Helper to extract a named numeric parameter from the proj4 string
    function get_param(name::String, default::Float64 = 0.0)
        m = match(Regex("\\+$(name)=([\\-\\+]?[\\d\\.eE\\+\\-]+)"), sr)
        return m === nothing ? default : parse(Float64, m.captures[1])
    end

    if proj == "longlat"
        return (
            GDTYP = Int32(1), P_ALP = 0.0, P_BET = 0.0, P_GAM = 0.0,
            XCENT = 0.0, YCENT = 0.0,
        )
    elseif proj == "lcc"
        return (
            GDTYP = Int32(2),
            P_ALP = get_param("lat_1"),
            P_BET = get_param("lat_2"),
            P_GAM = get_param("lon_0"),
            XCENT = get_param("lon_0"),
            YCENT = get_param("lat_0"),
        )
    elseif proj == "utm"
        return (
            GDTYP = Int32(5),
            P_ALP = get_param("zone"),
            P_BET = 0.0, P_GAM = 0.0, XCENT = 0.0, YCENT = 0.0,
        )
    elseif proj == "stere"
        return (
            GDTYP = Int32(6),
            P_ALP = get_param("lat_0"),
            P_BET = 0.0,
            P_GAM = get_param("lon_0"),
            XCENT = 0.0, YCENT = 0.0,
        )
    elseif proj == "merc"
        return (
            GDTYP = Int32(7),
            P_ALP = 0.0, P_BET = 0.0,
            P_GAM = get_param("lon_0"),
            XCENT = 0.0, YCENT = 0.0,
        )
    elseif proj == "tmerc"
        return (
            GDTYP = Int32(8),
            P_ALP = 0.0, P_BET = 0.0,
            P_GAM = get_param("lon_0"),
            XCENT = 0.0,
            YCENT = get_param("lat_0"),
        )
    else
        @warn "Unknown projection '$proj'; defaulting GDTYP to 1 (LATGRD3)"
        return (
            GDTYP = Int32(1), P_ALP = 0.0, P_BET = 0.0, P_GAM = 0.0,
            XCENT = 0.0, YCENT = 0.0,
        )
    end
end

"""
    _datetime_to_julian(dt::DateTime) -> Tuple{Int32, Int32}

Convert a `DateTime` to IOAPI Julian date and time.
Returns `(JDATE, JTIME)` where `JDATE = YYYYDDD` and `JTIME = HHMMSS`.
"""
function _datetime_to_julian(dt::DateTime)
    y = Dates.year(dt)
    doy = Dates.dayofyear(dt)
    jdate = Int32(y * 1000 + doy)
    jtime = Int32(Dates.hour(dt) * 10000 + Dates.minute(dt) * 100 + Dates.second(dt))
    return (jdate, jtime)
end

"""
    write_ioapi(filename::AbstractString, model_data::Dict{String, Array{Float64, 4}},
        grid::GridDef, hours::Vector{DateTime};
        n_layers::Int=1, description::String="", history::String="")

Write gridded emissions data to an IOAPI-format NetCDF file suitable for input to
air quality models (CMAQ, CAMx, etc.).

# Arguments
- `filename`: Output NetCDF file path.
- `model_data`: Dictionary mapping pollutant names to 4D arrays with dimensions
  `[nrows, ncols, n_layers, n_hours]`, as returned by `to_model_ready`.
- `grid::GridDef`: Grid definition (must be a regular grid with `Dx`, `Dy`, `X0`, `Y0`).
- `hours::Vector{DateTime}`: Ordered vector of time steps.
- `n_layers::Int=1`: Number of vertical layers.
- `description::String=""`: File description for the FILEDESC attribute.
- `history::String=""`: Processing history for the HISTORY attribute.

# IOAPI Format Details
The output file follows the IOAPI-3 conventions:
- Dimensions: `TSTEP` (unlimited), `LAY`, `ROW`, `COL`, `VAR`, `DATE-TIME` (=2)
- Coordinate variable `TFLAG(TSTEP, VAR, DATE-TIME)`: Julian date (YYYYDDD) and time (HHMMSS)
- One `Float32` variable per pollutant with dimensions `(COL, ROW, LAY, TSTEP)`
- Global attributes specify grid projection and structure (GDTYP, P_ALP, etc.)

# Example
```julia
grid = read_griddesc("GRIDDESC", "US36_CRO")
hours = DateTime(2019,7,1):Hour(1):DateTime(2019,7,1,23)
model_data = to_model_ready(merged, grid, collect(hours))
write_ioapi("emissions.nc", model_data, grid, collect(hours))
```
"""
function write_ioapi(
        filename::AbstractString, model_data::Dict{String, Array{Float64, 4}},
        grid::GridDef, hours::Vector{DateTime};
        n_layers::Int = 1, description::String = "", history::String = ""
    )
    if isempty(model_data)
        error("model_data is empty; nothing to write")
    end

    pollutant_names = sort(collect(keys(model_data)))
    nvars = length(pollutant_names)
    nrows = grid.Ny
    ncols = grid.Nx
    nsteps = length(hours)

    # Parse projection info
    proj_attrs = _parse_proj4_to_ioapi(grid.SR)

    # Compute TSTEP from hour spacing (default 1 hour)
    tstep = nsteps > 1 ? Int32(round(Dates.value(hours[2] - hours[1]) / 1000 / 3600) * 10000) : Int32(10000)

    # Start date/time from first hour
    sdate, stime = _datetime_to_julian(hours[1])

    # Build VAR-LIST: 16-char padded variable names
    varlist = join([rpad(v, 16) for v in pollutant_names])

    NCDatasets.Dataset(filename, "c") do ds
        # Define dimensions
        NCDatasets.defDim(ds, "TSTEP", Inf)  # unlimited
        NCDatasets.defDim(ds, "LAY", n_layers)
        NCDatasets.defDim(ds, "ROW", nrows)
        NCDatasets.defDim(ds, "COL", ncols)
        NCDatasets.defDim(ds, "VAR", nvars)
        NCDatasets.defDim(ds, "DATE-TIME", 2)

        # Global attributes
        ds.attrib["IOAPI_VERSION"] = "3.2"
        ds.attrib["EXEC_ID"] = "Emissions.jl"
        ds.attrib["FTYPE"] = Int32(1)  # GRDDED3
        ds.attrib["GDTYP"] = proj_attrs.GDTYP
        ds.attrib["P_ALP"] = Float64(proj_attrs.P_ALP)
        ds.attrib["P_BET"] = Float64(proj_attrs.P_BET)
        ds.attrib["P_GAM"] = Float64(proj_attrs.P_GAM)
        ds.attrib["XCENT"] = Float64(proj_attrs.XCENT)
        ds.attrib["YCENT"] = Float64(proj_attrs.YCENT)
        ds.attrib["XORIG"] = Float64(grid.X0)
        ds.attrib["YORIG"] = Float64(grid.Y0)
        ds.attrib["XCELL"] = Float64(grid.Dx)
        ds.attrib["YCELL"] = Float64(grid.Dy)
        ds.attrib["NCOLS"] = Int32(ncols)
        ds.attrib["NROWS"] = Int32(nrows)
        ds.attrib["NLAYS"] = Int32(n_layers)
        ds.attrib["NVARS"] = Int32(nvars)
        ds.attrib["SDATE"] = sdate
        ds.attrib["STIME"] = stime
        ds.attrib["TSTEP"] = tstep
        ds.attrib["NTHIK"] = Int32(1)
        ds.attrib["GDNAM"] = rpad(grid.Name, 16)
        ds.attrib["UPNAM"] = "Emissions.jl"
        ds.attrib["VAR-LIST"] = varlist
        ds.attrib["FILEDESC"] = description
        ds.attrib["HISTORY"] = history
        ds.attrib["VGTYP"] = Int32(-1)
        ds.attrib["VGTOP"] = Float32(0.0)
        ds.attrib["VGLVLS"] = Float32[1.0, 0.0]

        # TFLAG variable: (TSTEP, VAR, DATE-TIME)
        tflag_var = NCDatasets.defVar(
            ds, "TFLAG", Int32, ("DATE-TIME", "VAR", "TSTEP");
            attrib = Dict(
                "units" => "<YYYYDDD,HHMMSS>",
                "long_name" => rpad("TFLAG", 16),
                "var_desc" => rpad("Timestep-valid flags: <YYYYDDD,HHMMSS>", 80)
            )
        )

        # Write TFLAG data
        tflag_data = zeros(Int32, 2, nvars, nsteps)
        for t in 1:nsteps
            jdate, jtime = _datetime_to_julian(hours[t])
            for v in 1:nvars
                tflag_data[1, v, t] = jdate
                tflag_data[2, v, t] = jtime
            end
        end
        tflag_var[:, :, 1:nsteps] = tflag_data

        # Define and write pollutant variables: (COL, ROW, LAY, TSTEP)
        for (vi, polname) in enumerate(pollutant_names)
            pol_var = NCDatasets.defVar(
                ds, polname, Float32, ("COL", "ROW", "LAY", "TSTEP");
                attrib = Dict(
                    "units" => "mol/s or g/s",
                    "long_name" => rpad(polname, 16),
                    "var_desc" => rpad("Emissions of $polname", 80)
                )
            )

            arr = model_data[polname]
            # arr is [nrows, ncols, n_layers, n_hours]
            # NetCDF variable is (COL, ROW, LAY, TSTEP)
            out = zeros(Float32, ncols, nrows, n_layers, nsteps)
            for t in 1:nsteps, l in 1:n_layers, r in 1:nrows, c in 1:ncols
                out[c, r, l, t] = Float32(arr[r, c, l, t])
            end
            pol_var[:, :, :, 1:nsteps] = out
        end
    end

    return nothing
end

"""
    write_ioapi(filename::AbstractString, merged::DataFrame, grid::GridDef,
        hours::Vector{DateTime}; n_layers::Int=1, description::String="", history::String="")

Convenience method that accepts a merged `DataFrame` (from `merge_emissions`)
and converts it to model-ready arrays via `to_model_ready` before writing.

# Example
```julia
write_ioapi("emissions.nc", merged_df, grid, hours)
```
"""
function write_ioapi(
        filename::AbstractString, merged::DataFrame, grid::GridDef,
        hours::Vector{DateTime}; n_layers::Int = 1, description::String = "", history::String = ""
    )
    model_data = to_model_ready(merged, grid, hours; n_layers = n_layers)
    return write_ioapi(
        filename, model_data, grid, hours;
        n_layers = n_layers, description = description, history = history
    )
end
