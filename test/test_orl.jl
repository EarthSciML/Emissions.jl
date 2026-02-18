using DataFrames, Unitful

@testset "ORL format tests" begin
    @testset "ORLNonPointDataFrame" begin
        # Create synthetic ORL nonpoint data (21 columns)
        ncols = length(Emissions.ORL_NONPOINT_COLUMNS)
        df = DataFrame([Symbol("col$i") => [
            if i == 1; 1001            # FIPS
            elseif i == 2; 2101        # SCC
            elseif i == 7; "NOX"       # POLID
            elseif i == 8; 100.0       # ANN_VALUE (tons/yr)
            else; "test"
            end
        ] for i in 1:ncols])
        result = ORLNonPointDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
        @test result.df[1, :SCC] == "0000002101"
        @test result.df[1, :POLID] == "NOX"
        @test unit(result.df[1, :ANN_VALUE]) == u"kg/s"
        @test hasproperty(result.df, :COUNTRY)
    end

    @testset "ORLNonPointDataFrame wrong columns" begin
        df = DataFrame(a = [1], b = [2])
        @test_throws DimensionMismatch ORLNonPointDataFrame(df)
    end

    @testset "ORLPointDataFrame" begin
        ncols = length(Emissions.ORL_POINT_COLUMNS)
        df = DataFrame([Symbol("col$i") => [
            if i == 1; 1001            # FIPS
            elseif i == 7; 2101        # SCC
            elseif i == 10; 100.0      # STKHGT (ft)
            elseif i == 11; 10.0       # STKDIAM (ft)
            elseif i == 12; 500.0      # STKTEMP (°F)
            elseif i == 13; 1000.0     # STKFLOW (ft³/s)
            elseif i == 14; 50.0       # STKVEL (ft/s)
            elseif i == 19; -73.5      # LONGITUDE
            elseif i == 20; 40.5       # LATITUDE
            elseif i == 22; "NOX"      # POLID
            elseif i == 23; 100.0      # ANN_VALUE (tons/yr)
            else; "test"
            end
        ] for i in 1:ncols])
        result = ORLPointDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
        @test result.df[1, :SCC] == "0000002101"
        @test result.df[1, :STKHGT] isa Float64
        @test result.df[1, :STKHGT] ≈ 100.0 * 0.3048
        @test result.df[1, :STKDIAM] isa Float64
        @test result.df[1, :STKDIAM] ≈ 10.0 * 0.3048
        @test result.df[1, :STKTEMP] isa Float64
        @test result.df[1, :STKTEMP] ≈ (500.0 - 32.0) * 5.0 / 9.0 + 273.15
        @test result.df[1, :STKVEL] isa Float64
        @test result.df[1, :STKVEL] ≈ 50.0 * 0.3048
    end

    @testset "ORLPointDataFrame wrong columns" begin
        df = DataFrame(a = [1], b = [2])
        @test_throws DimensionMismatch ORLPointDataFrame(df)
    end

    @testset "ORLNonRoadDataFrame" begin
        ncols = length(Emissions.ORL_NONROAD_COLUMNS)
        df = DataFrame([Symbol("col$i") => [
            if i == 1; 1001
            elseif i == 2; 2101
            elseif i == 7; "NOX"
            elseif i == 8; 100.0
            else; "test"
            end
        ] for i in 1:ncols])
        result = ORLNonRoadDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
    end

    @testset "ORLOnRoadDataFrame" begin
        ncols = length(Emissions.ORL_ONROAD_COLUMNS)
        df = DataFrame([Symbol("col$i") => [
            if i == 1; 1001
            elseif i == 2; 2101
            elseif i == 7; "NOX"
            elseif i == 8; 100.0
            else; "test"
            end
        ] for i in 1:ncols])
        result = ORLOnRoadDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
    end

    @testset "ORLFireDataFrame" begin
        ncols = length(Emissions.ORL_FIRE_COLUMNS)
        df = DataFrame([Symbol("col$i") => [
            if i == 1; 1001            # FIPS
            elseif i == 4; 2101        # SCC
            elseif i == 7; 40.5        # LATITUDE
            elseif i == 8; -73.5       # LONGITUDE
            elseif i == 11; 8000.0     # HEATCONTENT
            elseif i == 12; "PM25"     # POLID
            elseif i == 13; 50.0       # ANN_VALUE
            else; "test"
            end
        ] for i in 1:ncols])
        result = ORLFireDataFrame(df)
        @test result.df[1, :FIPS] == "01001"
        @test result.df[1, :SCC] == "0000002101"
        @test result.df[1, :POLID] == "PM25"
        @test unit(result.df[1, :ANN_VALUE]) == u"kg/s"
    end

    @testset "ORLFireDataFrame wrong columns" begin
        df = DataFrame(a = [1], b = [2])
        @test_throws DimensionMismatch ORLFireDataFrame(df)
    end

    @testset "read_orl dispatch" begin
        # Create a temporary nonpoint ORL file
        ncols = length(Emissions.ORL_NONPOINT_COLUMNS)
        fields = String[]
        for i in 1:ncols
            if i == 1
                push!(fields, "1001")
            elseif i == 2
                push!(fields, "2101")
            elseif i == 7
                push!(fields, "NOX")
            elseif i == 8
                push!(fields, "100.0")
            else
                push!(fields, "test")
            end
        end
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# ORL header comment")
            println(io, join(fields, ","))
        end
        result = read_orl(tmpfile, :nonpoint)
        @test result isa ORLNonPointDataFrame
        @test nrow(result.df) == 1
        rm(tmpfile)
    end

    @testset "read_orl invalid format" begin
        @test_throws ArgumentError read_orl("dummy.csv", :invalid)
    end
end
