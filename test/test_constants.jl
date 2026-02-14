using Unitful

@testset "Constants tests" begin
    @testset "tonperyear" begin
        @test tonperyear ≈ 907.185u"kg" / 31_536_000u"s"
        @test unit(tonperyear) == u"kg/s"
    end

    @testset "tonpermonth" begin
        @test tonpermonth ≈ 907.185u"kg" / 2_628_288u"s"
        @test unit(tonpermonth) == u"kg/s"
    end

    @testset "foot" begin
        @test foot ≈ (1 / 3.28084)u"m"
        @test unit(foot) == u"m"
    end

    @testset "kelvin" begin
        # Freezing point: 32°F = 273.15 K
        @test kelvin(32.0) ≈ 273.15u"K"
        # Boiling point: 212°F = 373.15 K
        @test kelvin(212.0) ≈ 373.15u"K"
        # Absolute zero: -459.67°F = 0 K
        @test kelvin(-459.67) ≈ 0.0u"K" atol=0.01u"K"
    end

    @testset "Pollutants" begin
        @test Pollutants["VOC"] == "VOC"
        @test Pollutants["NOX"] == "NOX"
        @test Pollutants["NH3"] == "NH3"
        @test Pollutants["SO2"] == "SO2"
        @test Pollutants["PM25-PRI"] == "PM25"
        @test Pollutants["EVP__VOC"] == "VOC"
        @test Pollutants["EXH__PM25-PRI"] == "PM25"
        @test !haskey(Pollutants, "NONEXISTENT")
    end
end
