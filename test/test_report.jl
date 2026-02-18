using DataFrames, Dates, Unitful

@testset "Report tests" begin
    @testset "ReportConfig defaults" begin
        config = ReportConfig()
        @test config.group_by == [:POLID]
        @test config.time_resolution == :annual
        @test config.top_n == 0
    end

    @testset "ReportConfig invalid resolution" begin
        @test_throws ArgumentError ReportConfig(time_resolution = :invalid)
    end

    @testset "summary_by_pollutant pre-merge" begin
        df = DataFrame(
            FIPS = ["36001", "36001", "36005"],
            SCC = ["2103007000", "2103007000", "2103007000"],
            POLID = ["NOX", "NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0, 200.0]
        )
        result = summary_by_pollutant(df)
        @test nrow(result) == 2
        nox = filter(r -> r.pollutant == "NOX", result)
        voc = filter(r -> r.pollutant == "VOC", result)
        @test nox[1, :total] ≈ 150.0
        @test nox[1, :count] == 2
        @test voc[1, :total] ≈ 200.0
    end

    @testset "summary_by_pollutant post-merge" begin
        df = DataFrame(
            grid_row = [1, 1, 1],
            grid_col = [1, 1, 2],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "VOC", "NOX"],
            emission_rate = [100.0, 50.0, 25.0]
        )
        result = summary_by_pollutant(df)
        @test nrow(result) == 2
        nox = filter(r -> r.pollutant == "NOX", result)
        @test nox[1, :total] ≈ 125.0
    end

    @testset "summary_by_pollutant empty" begin
        df = DataFrame(POLID = String[], ANN_VALUE = Float64[])
        result = summary_by_pollutant(df)
        @test nrow(result) == 0
    end

    @testset "summary_by_pollutant with Unitful" begin
        df = DataFrame(
            POLID = ["NOX"],
            ANN_VALUE = [1.0e-3u"kg/s"]
        )
        result = summary_by_pollutant(df)
        @test nrow(result) == 1
        @test result[1, :total] ≈ 1.0e-3
    end

    @testset "summary_by_region" begin
        df = DataFrame(
            grid_row = [1, 1, 2],
            grid_col = [1, 2, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "NOX", "NOX"],
            emission_rate = [100.0, 50.0, 75.0]
        )
        region_map = Dict(
            (1, 1) => "East",
            (1, 2) => "East",
            (2, 1) => "West"
        )
        result = summary_by_region(df, region_map)
        @test nrow(result) == 2
        east = filter(r -> r.region == "East", result)
        west = filter(r -> r.region == "West", result)
        @test east[1, :total] ≈ 150.0
        @test west[1, :total] ≈ 75.0
    end

    @testset "summary_by_region empty" begin
        df = DataFrame(
            grid_row = Int[], grid_col = Int[],
            hour = DateTime[], pollutant = String[], emission_rate = Float64[]
        )
        result = summary_by_region(df, Dict{Tuple{Int, Int}, String}())
        @test nrow(result) == 0
    end

    @testset "summary_by_scc" begin
        df = DataFrame(
            FIPS = ["36001", "36001", "36005"],
            SCC = ["2103007000", "2103007000", "2103007001"],
            POLID = ["NOX", "NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0, 200.0]
        )
        result = summary_by_scc(df)
        @test nrow(result) == 2
        @test result[1, :total] ≈ 200.0  # sorted by total descending
    end

    @testset "summary_by_scc empty" begin
        df = DataFrame(SCC = String[], POLID = String[], ANN_VALUE = Float64[])
        result = summary_by_scc(df)
        @test nrow(result) == 0
    end

    @testset "summary_by_time daily" begin
        df = DataFrame(
            grid_row = [1, 1, 1],
            grid_col = [1, 1, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 12), DateTime(2019, 7, 2, 0)],
            pollutant = ["NOX", "NOX", "NOX"],
            emission_rate = [100.0, 50.0, 75.0]
        )
        result = summary_by_time(df; resolution = :daily)
        @test nrow(result) == 2
        day1 = filter(r -> r.period == "2019-07-01", result)
        @test day1[1, :total] ≈ 150.0
    end

    @testset "summary_by_time hourly" begin
        df = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)],
            pollutant = ["NOX", "NOX"],
            emission_rate = [100.0, 50.0]
        )
        result = summary_by_time(df; resolution = :hourly)
        @test nrow(result) == 2
    end

    @testset "summary_by_time monthly" begin
        df = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 8, 1, 0)],
            pollutant = ["NOX", "NOX"],
            emission_rate = [100.0, 50.0]
        )
        result = summary_by_time(df; resolution = :monthly)
        @test nrow(result) == 2
        jul = filter(r -> r.period == "2019-07", result)
        @test jul[1, :total] ≈ 100.0
    end

    @testset "summary_by_time invalid resolution" begin
        df = DataFrame(
            grid_row = [1], grid_col = [1],
            hour = [DateTime(2019, 7, 1)],
            pollutant = ["NOX"], emission_rate = [100.0]
        )
        @test_throws ArgumentError summary_by_time(df; resolution = :invalid)
    end

    @testset "compare_inventories" begin
        base = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007001"],
            POLID = ["NOX", "VOC"],
            ANN_VALUE = [100.0, 200.0]
        )
        scenario = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007001"],
            POLID = ["NOX", "VOC"],
            ANN_VALUE = [120.0, 180.0]
        )
        result = compare_inventories(base, scenario)
        @test nrow(result) == 2
        nox = filter(r -> r.POLID == "NOX", result)
        @test nox[1, :abs_diff] ≈ 20.0
        @test nox[1, :pct_diff] ≈ 20.0
        voc = filter(r -> r.POLID == "VOC", result)
        @test voc[1, :abs_diff] ≈ -20.0
        @test voc[1, :pct_diff] ≈ -10.0
    end

    @testset "compare_inventories new source" begin
        base = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0]
        )
        scenario = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007001"],
            POLID = ["NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0]
        )
        result = compare_inventories(base, scenario)
        @test nrow(result) == 2
        new_source = filter(r -> r.POLID == "VOC", result)
        @test new_source[1, :base_value] ≈ 0.0
        @test new_source[1, :scenario_value] ≈ 50.0
    end

    @testset "compare_inventories empty" begin
        base = DataFrame(FIPS = String[], SCC = String[], POLID = String[], ANN_VALUE = Float64[])
        scenario = DataFrame(FIPS = String[], SCC = String[], POLID = String[], ANN_VALUE = Float64[])
        result = compare_inventories(base, scenario)
        @test nrow(result) == 0
    end

    @testset "emissions_report basic" begin
        df = DataFrame(
            POLID = ["NOX", "NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0, 200.0]
        )
        result = emissions_report(df)
        @test nrow(result) == 2
        @test result[1, :total] ≈ 200.0  # VOC first (highest)
    end

    @testset "emissions_report with top_n" begin
        df = DataFrame(
            POLID = ["NOX", "VOC", "SO2"],
            ANN_VALUE = [100.0, 200.0, 50.0]
        )
        config = ReportConfig(top_n = 2)
        result = emissions_report(df; config = config)
        @test nrow(result) == 2
    end

    @testset "emissions_report post-merge" begin
        df = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "VOC"],
            emission_rate = [100.0, 50.0]
        )
        config = ReportConfig(group_by = [:pollutant])
        result = emissions_report(df; config = config)
        @test nrow(result) == 2
    end

    @testset "emissions_report empty" begin
        df = DataFrame()
        result = emissions_report(df)
        @test nrow(result) == 0
    end
end
