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
        @test unit(result.df[1, :STKHGT]) == u"m"
        @test unit(result.df[1, :STKTEMP]) == u"K"
        @test unit(result.df[1, :STKVEL]) == u"m/s"
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
