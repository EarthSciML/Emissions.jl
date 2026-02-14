export FF10NonPointDataFrame, FF10PointDataFrame, FF10NonRoadDataFrame, FF10OnRoadDataFrame

# Column names for FF10 nonpoint format.
# https://www.cmascenter.org/smoke/documentation/4.8.1/html/ch08s02s04.html#d0e38214
const FF10_NONPOINT_COLUMNS = [
    "COUNTRY", "FIPS", "TRIBAL_CODE", "CENSUS_TRACT", "SHAPE_ID", "SCC",
    "EMIS_TYPE", "POLID", "ANN_VALUE",
    "ANN_PCT_RED", "CONTROL_IDS", "CONTROL_MEASURES", "CURRENT_COST", "CUMULATIVE_COST", "PROJECTION_FACTOR",
    "REG_CODES", "CALC_METHOD", "CALC_YEAR", "DATE_UPDATED", "DATA_SET_ID", "JAN_VALUE", "FEB_VALUE", "MAR_VALUE",
    "APR_VALUE", "MAY_VALUE", "JUN_VALUE", "JUL_VALUE", "AUG_VALUE", "SEP_VALUE", "OCT_VALUE", "NOV_VALUE", "DEC_VALUE",
    "JAN_PCTRED", "FEB_PCTRED", "MAR_PCTRED", "APR_PCTRED", "MAY_PCTRED", "JUN_PCTRED", "JUL_PCTRED", "AUG_PCTRED",
    "SEP_PCTRED", "OCT_PCTRED", "NOV_PCTRED", "DEC_PCTRED", "COMMENT",
]

# Column names for FF10 point format.
# https://www.cmascenter.org/smoke/documentation/4.8.1/html/ch08s02s08.html#sect_input_ptinv_ff10
const FF10_POINT_COLUMNS = [
    "COUNTRY", "FIPS", "TRIBAL_CODE", "FACILITY_ID", "UNIT_ID", "REL_POINT_ID", "PROCESS_ID", "AGY_FACILITY_ID", "AGY_UNIT_ID",
    "AGY_REL_POINT_ID", "AGY_PROCESS_ID", "SCC", "POLID", "ANN_VALUE", "ANN_PCT_RED", "FACILITY_NAME", "ERPTYPE", "STKHGT",
    "STKDIAM", "STKTEMP", "STKFLOW", "STKVEL", "NAICS", "LONGITUDE", "LATITUDE", "LL_DATUM", "HORIZ_COLL_MTHD", "DESIGN_CAPACITY",
    "DESIGN_CAPACITY_UNITS", "REG_CODES", "FAC_SOURCE_TYPE", "UNIT_TYPE_CODE", "CONTROL_IDS", "CONTROL_MEASURES",
    "CURRENT_COST", "CUMULATIVE_COST", "PROJECTION_FACTOR", "SUBMITTER_FAC_ID", "CALC_METHOD", "DATA_SET_ID", "FACIL_CATEGORY_CODE",
    "ORIS_FACILITY_CODE", "ORIS_BOILER_ID", "IPM_YN", "CALC_YEAR", "DATE_UPDATED", "FUG_HEIGHT", "FUG_WIDTH_XDIM", "FUG_LENGTH_YDIM",
    "FUG_ANGLE", "ZIPCODE", "ANNUAL_AVG_HOURS_PER_YEAR", "JAN_VALUE", "FEB_VALUE", "MAR_VALUE",
    "APR_VALUE", "MAY_VALUE", "JUN_VALUE", "JUL_VALUE", "AUG_VALUE", "SEP_VALUE", "OCT_VALUE", "NOV_VALUE", "DEC_VALUE",
    "JAN_PCTRED", "FEB_PCTRED", "MAR_PCTRED", "APR_PCTRED", "MAY_PCTRED", "JUN_PCTRED", "JUL_PCTRED", "AUG_PCTRED",
    "SEP_PCTRED", "OCT_PCTRED", "NOV_PCTRED", "DEC_PCTRED", "COMMENT",
]

const MONTHLY_VALUE_COLUMNS = [
    "JAN_VALUE", "FEB_VALUE", "MAR_VALUE",
    "APR_VALUE", "MAY_VALUE", "JUN_VALUE",
    "JUL_VALUE", "AUG_VALUE", "SEP_VALUE",
    "OCT_VALUE", "NOV_VALUE", "DEC_VALUE",
]

"""
    transform_fips!(df::DataFrame)

Standardize FIPS codes in a DataFrame: strip leading country digit if present,
then left-pad to 5 characters.
"""
function transform_fips!(df::DataFrame)
    fips_transformed = String[]
    for fips in df[!, :FIPS]
        fips_str = string(fips)
        if length(fips_str) == 6
            fips_str = fips_str[2:end]
        end
        push!(fips_transformed, lpad(fips_str, 5, '0'))
    end
    df[!, :FIPS] = fips_transformed
    return df
end

"""
    convert_emissions_units!(df::DataFrame)

Convert annual and monthly emission values from tons to kg/s.
"""
function convert_emissions_units!(df::DataFrame)
    df.ANN_VALUE = df.ANN_VALUE * tonperyear
    for col in MONTHLY_VALUE_COLUMNS
        df[!, col] = df[!, col] * tonpermonth
    end
    return df
end

"""
    FF10NonPointDataFrame <: EmissionsDataFrame

Wrapper for EPA FF10 nonpoint emissions data. Validates that the DataFrame has 45 columns,
renames them to standard names, standardizes FIPS/SCC codes, and converts units.
"""
struct FF10NonPointDataFrame <: EmissionsDataFrame
    df::DataFrame

    FF10NonPointDataFrame(df::DataFrame) = begin
        if size(df, 2) != 45
            throw(DimensionMismatch("FF10 nonpoint file should have 45 fields but instead has $(size(df,2))"))
        end

        rename!(df, FF10_NONPOINT_COLUMNS)
        transform_fips!(df)
        df[!, :SCC] = [lpad(scc, 10, '0') for scc in string.(df[!, :SCC])]
        convert_emissions_units!(df)

        return new(df)
    end
end

"""
    FF10NonRoadDataFrame <: EmissionsDataFrame

Wrapper for EPA FF10 nonroad emissions data. Uses the same format as nonpoint.
"""
struct FF10NonRoadDataFrame <: EmissionsDataFrame
    df::DataFrame

    FF10NonRoadDataFrame(df::DataFrame) = new(FF10NonPointDataFrame(df).df)
end

"""
    FF10OnRoadDataFrame <: EmissionsDataFrame

Wrapper for EPA FF10 on-road emissions data. Uses the same format as nonpoint.
"""
struct FF10OnRoadDataFrame <: EmissionsDataFrame
    df::DataFrame

    FF10OnRoadDataFrame(df::DataFrame) = new(FF10NonPointDataFrame(df).df)
end

"""
    FF10PointDataFrame <: EmissionsDataFrame

Wrapper for EPA FF10 point source emissions data. Validates that the DataFrame has 77 columns,
renames them to standard names, standardizes FIPS/SCC codes, converts units,
and converts stack parameters from imperial to SI units.
"""
struct FF10PointDataFrame <: EmissionsDataFrame
    df::DataFrame

    FF10PointDataFrame(df::DataFrame) = begin
        if size(df, 2) != 77
            throw(DimensionMismatch("FF10 point file should have 77 fields but instead has $(size(df,2))"))
        end

        rename!(df, FF10_POINT_COLUMNS)
        transform_fips!(df)
        df[!, :SCC] = [lpad(scc, 10, '0') for scc in string.(df[!, :SCC])]
        convert_emissions_units!(df)

        # Convert stack parameters from imperial to SI
        df.STKHGT = df.STKHGT * foot
        df.STKDIAM = df.STKDIAM * foot
        df.STKTEMP = kelvin.(df.STKTEMP)
        df.STKFLOW = df.STKFLOW * foot * foot * foot / u"s"
        df.STKVEL = df.STKVEL * foot / u"s"

        return new(df)
    end
end
