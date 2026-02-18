export writeEmis, find_surrogate_by_code, get_data_weight_shapefiles

"""
    format_float(val::Number)

Format a floating point number for output. Uses fixed-point notation
for small numbers and scientific notation for large ones.
"""
function format_float(val::Number)
    if abs(val) < 1e-20
        return "0.0"
    elseif abs(val) < 0.001 || abs(val) > 1e6
        return @sprintf("%.6e", val)
    else
        return @sprintf("%.6f", val)
    end
end

"""
    find_surrogate_by_code(srgSpecs::Vector{SurrogateSpec}, code::Int; region::String="") -> Union{SurrogateSpec, Nothing}

Find a surrogate specification by its code number, with optional region filtering.

When `region` is non-empty, only matches surrogates with that region.
When `region` is empty (default), matches any region.

Returns the matching `SurrogateSpec` or `nothing` if not found.
"""
function find_surrogate_by_code(srgSpecs::Vector{SurrogateSpec}, code::Int; region::String = "")
    for srg in srgSpecs
        if srg.Code == code && (isempty(region) || srg.Region == region)
            return srg
        end
    end
    return nothing
end

# Keep the two-argument region version for backward compatibility
"""
    find_surrogate_by_code(srgSpecs::Vector{SurrogateSpec}, region::String, code::Int)

Find a surrogate specification by its region and code number.
Backward-compatible method; prefer using the keyword argument version.
"""
function find_surrogate_by_code(srgSpecs::Vector{SurrogateSpec}, region::String, code::Int)
    return find_surrogate_by_code(srgSpecs, code; region = region)
end

"""
    get_data_weight_shapefiles(srg::SurrogateSpec)

Return a tuple of (data_shapefile, weight_shapefile) paths from a surrogate specification.
"""
function get_data_weight_shapefiles(srg::SurrogateSpec)
    return (srg.DataShapefile, srg.WeightShapefile)
end

"""
    writeEmis(filename, grid_data, grid::GridDef; pollutant="", units="")

Write gridded emissions data to a CSV file.

# Arguments
- `filename`: Output file path
- `grid_data`: Matrix or sparse matrix of emissions values
- `grid`: Grid definition
- `pollutant`: Pollutant name (optional, for header)
- `units`: Units string (optional, for header)
"""
function writeEmis(filename::AbstractString, grid_data, grid::GridDef;
    pollutant::AbstractString="", units::AbstractString="")

    open(filename, "w") do io
        # Write header
        println(io, "# Gridded emissions: $pollutant ($units)")
        println(io, "# Grid: $(grid.Name) $(grid.Nx)x$(grid.Ny)")
        println(io, "row,col,value")

        for j in 1:grid.Ny
            for i in 1:grid.Nx
                val = grid_data[j, i]
                if val != 0.0
                    println(io, "$j,$i,$(format_float(val))")
                end
            end
        end
    end
end
