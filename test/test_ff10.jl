using DataFrames, Unitful

@testset "FF10 tests" begin
    @testset "FF10NonPointDataFrame" begin
        # Create a synthetic 45-column DataFrame
        df = DataFrame([Symbol("col$i") => [i == 9 ? 100.0 : (i in 21:32 ? 10.0 : (i == 2 ? 1001 : (i == 6 ? 2101 : "test")))] for i in 1:45])
        result = FF10NonPointDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
        @test result.df[1, :SCC] == "0000002101"
        @test unit(result.df[1, :ANN_VALUE]) == u"kg/s"
    end

    @testset "FF10NonPointDataFrame wrong columns" begin
        df = DataFrame(a=[1], b=[2])
        @test_throws DimensionMismatch FF10NonPointDataFrame(df)
    end

    @testset "FF10PointDataFrame" begin
        # Create a synthetic 77-column DataFrame
        df = DataFrame([Symbol("col$i") => [
            if i == 2; 1001          # FIPS
            elseif i == 12; 2101     # SCC
            elseif i == 14; 100.0    # ANN_VALUE
            elseif i == 18; 100.0    # STKHGT (ft)
            elseif i == 19; 10.0     # STKDIAM (ft)
            elseif i == 20; 500.0    # STKTEMP (°F)
            elseif i == 21; 1000.0   # STKFLOW (ft³/s)
            elseif i == 22; 50.0     # STKVEL (ft/s)
            elseif i in 53:64; 10.0  # Monthly values
            else; "test"
            end
        ] for i in 1:77])
        result = FF10PointDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
        # Stack params should be plain Float64 (no Unitful), ready for downstream use
        @test result.df[1, :STKHGT] isa Float64
        @test result.df[1, :STKHGT] ≈ 100.0 * 0.3048  # 100 ft in meters
        @test result.df[1, :STKDIAM] isa Float64
        @test result.df[1, :STKDIAM] ≈ 10.0 * 0.3048
        @test result.df[1, :STKTEMP] isa Float64
        @test result.df[1, :STKTEMP] ≈ (500.0 - 32.0) * 5.0 / 9.0 + 273.15  # 500°F in K
        @test result.df[1, :STKVEL] isa Float64
        @test result.df[1, :STKVEL] ≈ 50.0 * 0.3048
    end

    @testset "FF10 point source end-to-end with classify and laypoint" begin
        # Create point data with realistic stack parameters
        df = DataFrame([Symbol("col$i") => [
            if i == 2; 1001              # FIPS
            elseif i == 12; 2101         # SCC
            elseif i == 13; "NOX"        # POLID
            elseif i == 14; 100.0        # ANN_VALUE
            elseif i == 18; 200.0        # STKHGT (ft) = 60.96m
            elseif i == 19; 10.0         # STKDIAM (ft)
            elseif i == 20; 500.0        # STKTEMP (°F)
            elseif i == 21; 1000.0       # STKFLOW (ft³/s)
            elseif i == 22; 50.0         # STKVEL (ft/s)
            elseif i == 24; -73.5        # LONGITUDE
            elseif i == 25; 40.5         # LATITUDE
            elseif i in 53:64; 10.0      # Monthly values
            else; "test"
            end
        ] for i in 1:77])
        point_df = FF10PointDataFrame(df)

        # classify_point_sources should work with plain Float64 stack params
        classified = classify_point_sources(point_df.df)
        @test classified[1, :source_class] in ["surface", "elevated", "ping"]
        @test classified[1, :analytical_plume_rise] isa Float64
        @test classified[1, :analytical_plume_rise] >= 0.0

        # laypoint should work with plain Float64 stack params
        layer_config = LayerConfig([0.0, 50.0, 100.0, 200.0, 500.0, 1000.0])
        met = MetProfile(
            fill(293.15, 5),   # temperature
            fill(5.0, 5),      # wind_speed
            fill(0.25, 5),     # stability_class (unstable)
            fill(0.01, 5),     # stability_param
        )
        met_profiles = Dict("default" => met)
        result = laypoint(point_df.df, met_profiles, layer_config)
        @test nrow(result) >= 1
        @test hasproperty(result, :layer)
        @test hasproperty(result, :layer_fraction)
        @test sum(result.layer_fraction) ≈ 1.0
    end

    @testset "FF10PointDataFrame wrong columns" begin
        df = DataFrame(a=[1], b=[2])
        @test_throws DimensionMismatch FF10PointDataFrame(df)
    end

    @testset "FF10NonRoadDataFrame" begin
        df = DataFrame([Symbol("col$i") => [i == 9 ? 100.0 : (i in 21:32 ? 10.0 : (i == 2 ? 1001 : (i == 6 ? 2101 : "test")))] for i in 1:45])
        result = FF10NonRoadDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
    end

    @testset "FF10OnRoadDataFrame" begin
        df = DataFrame([Symbol("col$i") => [i == 9 ? 100.0 : (i in 21:32 ? 10.0 : (i == 2 ? 1001 : (i == 6 ? 2101 : "test")))] for i in 1:45])
        result = FF10OnRoadDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
    end

    @testset "transform_fips!" begin
        df = DataFrame(FIPS=[100101, 1001, 1])
        Emissions.transform_fips!(df)
        @test df.FIPS[1] == "00101"
        @test df.FIPS[2] == "01001"
        @test df.FIPS[3] == "00001"
    end

    @testset "FIPS with 6 digits" begin
        df = DataFrame(FIPS=[901001])
        Emissions.transform_fips!(df)
        @test df.FIPS[1] == "01001"
    end
end
