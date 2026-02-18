export ORLNonPointDataFrame, ORLPointDataFrame, ORLNonRoadDataFrame,
    ORLOnRoadDataFrame, ORLFireDataFrame

# Column names for ORL nonpoint format.
# https://www.cmascenter.org/smoke/documentation/4.8.1/html/ch08s02s03.html
const ORL_NONPOINT_COLUMNS = [
    "FIPS", "SCC", "SIC", "MACT", "SRCTYPE", "NAICS", "POLID",
    "ANN_VALUE", "AVG_DAY_VALUE", "CONTROL_EFF", "RULE_EFF",
    "RULE_PEN", "EMIS_TYPE", "CALC_METHOD", "CALC_YEAR",
    "DATE_UPDATED", "DATA_SET_ID", "CEFF", "REFF", "RPEN", "COMMENT",
]

# Column names for ORL point format.
const ORL_POINT_COLUMNS = [
    "FIPS", "PLANTID", "POINTID", "STACKID", "SEGMENT", "PLANT",
    "SCC", "ERPTYPE", "SRCTYPE", "STKHGT", "STKDIAM", "STKTEMP",
    "STKFLOW", "STKVEL", "SIC", "MACT", "NAICS", "CTYPE", "LONGITUDE",
    "LATITUDE", "UTMZ", "POLID", "ANN_VALUE", "AVG_DAY_VALUE",
    "CONTROL_EFF", "RULE_EFF", "RULE_PEN", "EMIS_TYPE",
    "CALC_METHOD", "CALC_YEAR", "DATE_UPDATED", "DATA_SET_ID",
    "CEFF", "REFF", "RPEN", "CPRI", "CSEC", "COMMENT",
]

# Column names for ORL nonroad format (same as nonpoint).
const ORL_NONROAD_COLUMNS = ORL_NONPOINT_COLUMNS

# Column names for ORL onroad format (same as nonpoint).
const ORL_ONROAD_COLUMNS = ORL_NONPOINT_COLUMNS

# Column names for ORL fire format (PTFIRE).
const ORL_FIRE_COLUMNS = [
    "FIPS", "FIREID", "LOCID", "SCC", "DATA_SOURCE", "FIRENAME",
    "LATITUDE", "LONGITUDE", "NFDRSCODE", "MATBURNED", "HEATCONTENT",
    "POLID", "ANN_VALUE", "AVG_DAY_VALUE", "DATE_START", "DATE_END",
    "FACILID", "AESSION_NBR", "DATE_UPDATED", "DATA_SET_ID", "COMMENT",
]

"""
    ORLNonPointDataFrame <: EmissionsDataFrame

Wrapper for legacy ORL nonpoint emissions data. Validates that the DataFrame has
$(length(ORL_NONPOINT_COLUMNS)) columns, renames them to standard names, standardizes
FIPS/SCC codes, and converts units from tons/year to kg/s.
"""
struct ORLNonPointDataFrame <: EmissionsDataFrame
    df::DataFrame

    ORLNonPointDataFrame(df::DataFrame) = begin
        ncols = length(ORL_NONPOINT_COLUMNS)
        if size(df, 2) != ncols
            throw(DimensionMismatch(
                "ORL nonpoint file should have $ncols fields but instead has $(size(df,2))"
            ))
        end

        rename!(df, ORL_NONPOINT_COLUMNS)
        transform_fips!(df)
        df[!, :SCC] = [lpad(scc, 10, '0') for scc in string.(df[!, :SCC])]
        # Convert from tons/year to kg/s
        df.ANN_VALUE = df.ANN_VALUE * tonperyear

        # Add COUNTRY column (derived from original FIPS)
        if !hasproperty(df, :COUNTRY)
            df[!, :COUNTRY] .= "US"
        end

        return new(df)
    end
end

"""
    ORLNonRoadDataFrame <: EmissionsDataFrame

Wrapper for legacy ORL nonroad emissions data. Uses the same format as ORL nonpoint.
"""
struct ORLNonRoadDataFrame <: EmissionsDataFrame
    df::DataFrame

    ORLNonRoadDataFrame(df::DataFrame) = new(ORLNonPointDataFrame(df).df)
end

"""
    ORLOnRoadDataFrame <: EmissionsDataFrame

Wrapper for legacy ORL on-road emissions data. Uses the same format as ORL nonpoint.
"""
struct ORLOnRoadDataFrame <: EmissionsDataFrame
    df::DataFrame

    ORLOnRoadDataFrame(df::DataFrame) = new(ORLNonPointDataFrame(df).df)
end

"""
    ORLPointDataFrame <: EmissionsDataFrame

Wrapper for legacy ORL point source emissions data. Validates that the DataFrame has
$(length(ORL_POINT_COLUMNS)) columns, renames them to standard names, standardizes
FIPS/SCC codes, converts emission units, and converts stack parameters from imperial
to SI units (feet→m, °F→K).
"""
struct ORLPointDataFrame <: EmissionsDataFrame
    df::DataFrame

    ORLPointDataFrame(df::DataFrame) = begin
        ncols = length(ORL_POINT_COLUMNS)
        if size(df, 2) != ncols
            throw(DimensionMismatch(
                "ORL point file should have $ncols fields but instead has $(size(df,2))"
            ))
        end

        rename!(df, ORL_POINT_COLUMNS)
        transform_fips!(df)
        df[!, :SCC] = [lpad(scc, 10, '0') for scc in string.(df[!, :SCC])]
        df.ANN_VALUE = df.ANN_VALUE * tonperyear

        # Convert stack parameters from imperial to SI
        df.STKHGT = ustrip.(u"m", df.STKHGT * foot)
        df.STKDIAM = ustrip.(u"m", df.STKDIAM * foot)
        df.STKTEMP = ustrip.(u"K", kelvin.(df.STKTEMP))
        df.STKFLOW = ustrip.(u"m^3/s", df.STKFLOW * foot^3 / u"s")
        df.STKVEL = ustrip.(u"m/s", df.STKVEL * foot / u"s")

        if !hasproperty(df, :COUNTRY)
            df[!, :COUNTRY] .= "US"
        end

        return new(df)
    end
end

"""
    ORLFireDataFrame <: EmissionsDataFrame

Wrapper for legacy ORL fire emissions data (PTFIRE format). Validates that the DataFrame has
$(length(ORL_FIRE_COLUMNS)) columns, renames them to standard names, standardizes
FIPS/SCC codes, and converts emission units.

Includes fire-specific fields: LATITUDE, LONGITUDE, HEATCONTENT, DATE_START, DATE_END.
"""
struct ORLFireDataFrame <: EmissionsDataFrame
    df::DataFrame

    ORLFireDataFrame(df::DataFrame) = begin
        ncols = length(ORL_FIRE_COLUMNS)
        if size(df, 2) != ncols
            throw(DimensionMismatch(
                "ORL fire file should have $ncols fields but instead has $(size(df,2))"
            ))
        end

        rename!(df, ORL_FIRE_COLUMNS)
        transform_fips!(df)
        df[!, :SCC] = [lpad(scc, 10, '0') for scc in string.(df[!, :SCC])]
        df.ANN_VALUE = df.ANN_VALUE * tonperyear

        if !hasproperty(df, :COUNTRY)
            df[!, :COUNTRY] .= "US"
        end

        return new(df)
    end
end
