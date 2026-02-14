export tonperyear, tonpermonth, foot, kelvin, Pollutants

const tonperyear = 907.185u"kg" / 31_536_000u"s"
const tonpermonth = 907.185u"kg" / 2_628_288u"s"
const foot = (1 / 3.28084)u"m"

"""
    kelvin(F)

Convert temperature from Fahrenheit to Kelvin.
"""
kelvin(F) = ((F - 32.0) * 5.0 / 9.0 + 273.15) * u"K"

"""
    Pollutants

Dictionary mapping specific pollutant identifiers to aggregated pollutant group names.
"""
const Pollutants = Dict(
    "EVP__VOC" => "VOC", "EXH__VOC" => "VOC", "VOC" => "VOC", "VOC_INV" => "VOC",
    "NOX" => "NOX",
    "NH3" => "NH3",
    "SO2" => "SO2",
    "BRK__PM25-PRI" => "PM25", "EXH__PM25-PRI" => "PM25", "EXH__PM2_5" => "PM25",
    "PM25-PRI" => "PM25", "PM25TOTAL" => "PM25", "PM2_5" => "PM25", "TIR__PM25-PRI" => "PM25",
    "BRK__PM2_5" => "PM25", "TIR__PM2_5" => "PM25",
)
