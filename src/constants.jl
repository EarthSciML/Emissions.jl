export tonperyear, tonpermonth, foot, kelvin, Pollutants

"""
    tonperyear

Conversion factor from tons per year to kg per second.
1 ton = 907.185 kg, 1 year = 31,536,000 seconds.
"""
const tonperyear = 907.185u"kg" / 31_536_000u"s"

"""
    tonpermonth

Conversion factor from tons per month to kg per second.
1 ton = 907.185 kg, 1 month = 2,628,288 seconds (average month).
"""
const tonpermonth = 907.185u"kg" / 2_628_288u"s"

"""
    foot

Conversion factor from feet to meters.
1 foot = 0.3048 meters.
"""
const foot = 0.3048u"m"

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
