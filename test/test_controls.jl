using DataFrames

@testset "Controls tests" begin
    @testset "read_growth_factors" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Growth factors file")
            println(io, "036001;2103007000;VOC;2017;2025;1.5")
            println(io, "000000;2103007000;NOX;2017;2025;0.8 !national default")
        end
        controls = read_growth_factors(tmpfile)
        @test length(controls) == 2
        @test controls[1].fips == "36001"
        @test controls[1].scc == "2103007000"
        @test controls[1].pollutant == "VOC"
        @test controls[1].control_type == :growth
        @test controls[1].factor ≈ 1.5
        @test controls[1].base_year == 2017
        @test controls[1].target_year == 2025
        @test controls[2].fips == "00000"
        @test controls[2].factor ≈ 0.8
        rm(tmpfile)
    end

    @testset "read_growth_factors empty file" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Only comments")
        end
        controls = read_growth_factors(tmpfile)
        @test isempty(controls)
        rm(tmpfile)
    end

    @testset "read_control_factors" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Control factors file")
            # efficiency=80%, effectiveness=90%, penetration=100%
            # factor = 1 - (80/100 * 90/100 * 100/100) = 1 - 0.72 = 0.28
            println(io, "036001;2103007000;VOC;80;90;100")
            # efficiency=50%, effectiveness=100%, penetration=50%
            # factor = 1 - (50/100 * 100/100 * 50/100) = 1 - 0.25 = 0.75
            println(io, "000000;0000000000;NOX;50;100;50")
        end
        controls = read_control_factors(tmpfile)
        @test length(controls) == 2
        @test controls[1].fips == "36001"
        @test controls[1].control_type == :multiplicative
        @test controls[1].factor ≈ 0.28
        @test controls[2].fips == "00000"
        @test controls[2].scc == "0000000000"
        @test controls[2].factor ≈ 0.75
        rm(tmpfile)
    end

    @testset "_match_control hierarchical matching" begin
        controls = [
            ControlSpec("", "36001", "2103007000", "VOC", :multiplicative, 0.5, 0, 0),  # Level 1
            ControlSpec("", "36001", "2103007000", "", :multiplicative, 0.6, 0, 0),       # Level 2
            ControlSpec("", "36001", "0000000000", "NOX", :multiplicative, 0.7, 0, 0),    # Level 3
            ControlSpec("", "36001", "0000000000", "", :multiplicative, 0.8, 0, 0),       # Level 4
            ControlSpec("", "00000", "2103007000", "SO2", :multiplicative, 0.3, 0, 0),    # Level 5
            ControlSpec("", "00000", "0000000000", "PM25", :multiplicative, 0.4, 0, 0),   # Level 6
            ControlSpec("", "00000", "0000000000", "", :multiplicative, 0.9, 0, 0),       # Level 7
        ]

        # Level 1: Exact FIPS + SCC + pollutant
        result = Emissions._match_control(controls, "36001", "2103007000", "VOC", :multiplicative)
        @test result !== nothing
        @test result.factor ≈ 0.5

        # Level 2: FIPS + SCC, any pollutant
        result = Emissions._match_control(controls, "36001", "2103007000", "CO", :multiplicative)
        @test result !== nothing
        @test result.factor ≈ 0.6

        # Level 3: FIPS + pollutant, any SCC
        result = Emissions._match_control(controls, "36001", "9999999999", "NOX", :multiplicative)
        @test result !== nothing
        @test result.factor ≈ 0.7

        # Level 4: FIPS only
        result = Emissions._match_control(controls, "36001", "9999999999", "CO", :multiplicative)
        @test result !== nothing
        @test result.factor ≈ 0.8

        # Level 5: National + SCC + pollutant
        result = Emissions._match_control(controls, "99999", "2103007000", "SO2", :multiplicative)
        @test result !== nothing
        @test result.factor ≈ 0.3

        # Level 6: National + pollutant
        result = Emissions._match_control(controls, "99999", "9999999999", "PM25", :multiplicative)
        @test result !== nothing
        @test result.factor ≈ 0.4

        # Level 7: National default
        result = Emissions._match_control(controls, "99999", "9999999999", "CO", :multiplicative)
        @test result !== nothing
        @test result.factor ≈ 0.9

        # No match for different control type
        result = Emissions._match_control(controls, "36001", "2103007000", "VOC", :growth)
        @test result === nothing
    end

    @testset "apply_controls growth only" begin
        emissions = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["VOC", "VOC"],
            ANN_VALUE = [100.0, 200.0],
        )
        controls = [
            ControlSpec("", "00000", "2103007000", "VOC", :growth, 1.5, 2017, 2025),
        ]

        result = apply_controls(emissions, controls)
        @test result[1, :ANN_VALUE] ≈ 150.0
        @test result[2, :ANN_VALUE] ≈ 300.0
        # Original should be unchanged
        @test emissions[1, :ANN_VALUE] ≈ 100.0
    end

    @testset "apply_controls multiplicative only" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0],
        )
        controls = [
            ControlSpec("", "00000", "2103007000", "NOX", :multiplicative, 0.5, 0, 0),
        ]

        result = apply_controls(emissions, controls)
        @test result[1, :ANN_VALUE] ≈ 50.0
    end

    @testset "apply_controls combined growth and multiplicative" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["VOC"],
            ANN_VALUE = [100.0],
        )
        controls = [
            ControlSpec("", "00000", "2103007000", "VOC", :growth, 2.0, 2017, 2025),
            ControlSpec("", "00000", "2103007000", "VOC", :multiplicative, 0.5, 0, 0),
        ]

        result = apply_controls(emissions, controls)
        # growth * multiplicative = 2.0 * 0.5 = 1.0
        @test result[1, :ANN_VALUE] ≈ 100.0
    end

    @testset "apply_controls no matching control" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["9999999999"],
            POLID = ["CO"],
            ANN_VALUE = [100.0],
        )
        controls = ControlSpec[]

        result = apply_controls(emissions, controls)
        @test result[1, :ANN_VALUE] ≈ 100.0
    end

    @testset "apply_controls zero factor" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["VOC"],
            ANN_VALUE = [100.0],
        )
        controls = [
            ControlSpec("", "00000", "2103007000", "VOC", :multiplicative, 0.0, 0, 0),
        ]

        result = apply_controls(emissions, controls)
        @test result[1, :ANN_VALUE] ≈ 0.0
    end

    @testset "apply_controls factor greater than 1" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["VOC"],
            ANN_VALUE = [100.0],
        )
        controls = [
            ControlSpec("", "00000", "2103007000", "VOC", :growth, 3.0, 2017, 2030),
        ]

        result = apply_controls(emissions, controls)
        @test result[1, :ANN_VALUE] ≈ 300.0
    end

    @testset "apply_controls preserves other columns" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["VOC"],
            ANN_VALUE = [100.0],
            COUNTRY = ["USA"],
            LONGITUDE = [-73.5],
        )
        controls = [
            ControlSpec("", "00000", "2103007000", "VOC", :growth, 2.0, 2017, 2025),
        ]

        result = apply_controls(emissions, controls)
        @test hasproperty(result, :COUNTRY)
        @test result[1, :COUNTRY] == "USA"
        @test result[1, :LONGITUDE] ≈ -73.5
        @test result[1, :ANN_VALUE] ≈ 200.0
    end
end
