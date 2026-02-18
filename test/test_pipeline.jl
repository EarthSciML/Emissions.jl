using DataFrames, CSV, Unitful, Emissions

@testset "Pipeline tests" begin
    @testset "normalize_country" begin
        @test normalize_country("US") == "USA"
        @test normalize_country("0") == "USA"
        @test normalize_country("1") == "Canada"
        @test normalize_country("2") == "Mexico"
        @test normalize_country("USA") == "USA"
        @test normalize_country("Canada") == "Canada"
        @test normalize_country(" US ") == "USA"
    end

    @testset "read_gridref" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# This is a comment")
            println(io, "036001;2103007000;100!population surrogate")
            println(io, "036005;2103007000;200")
            println(io, "000000;0000000000;100!fallback entry")
        end
        result = read_gridref(tmpfile)
        @test size(result, 1) == 3
        @test result.COUNTRY[1] == "USA"
        @test result.FIPS[1] == "36001"
        @test result.SCC[1] == "2103007000"
        @test result.Surrogate[1] == 100
        @test result.COUNTRY[2] == "USA"
        @test result.FIPS[2] == "36005"
        @test result.Surrogate[2] == 200
        # Fallback entry with FIPS "00000"
        @test result.FIPS[3] == "00000"
        @test result.Surrogate[3] == 100
        rm(tmpfile)
    end

    @testset "read_gridref skips comments and empty lines" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Header")
            println(io, "")
            println(io, "# Another comment")
            println(io, "036001;2103007000;100")
        end
        result = read_gridref(tmpfile)
        @test size(result, 1) == 1
        rm(tmpfile)
    end

    @testset "read_ff10 nonpoint" begin
        # Create a temp CSV with 45 columns
        tmpfile = tempname() * ".csv"
        row_data = []
        for i in 1:45
            if i == 2
                push!(row_data, "36001")  # FIPS
            elseif i == 6
                push!(row_data, "2103007000")  # SCC
            elseif i == 9
                push!(row_data, "100.0")  # ANN_VALUE
            elseif i in 21:32
                push!(row_data, "10.0")  # Monthly values
            else
                push!(row_data, "test")
            end
        end
        open(tmpfile, "w") do io
            println(io, join(row_data, ","))
        end
        result = read_ff10(tmpfile, :nonpoint)
        @test result isa FF10NonPointDataFrame
        @test size(result.df, 1) == 1
        @test result.df[1, :FIPS] == "36001"
        rm(tmpfile)
    end

    @testset "read_ff10 point" begin
        tmpfile = tempname() * ".csv"
        row_data = []
        for i in 1:77
            if i == 2
                push!(row_data, "36001")  # FIPS
            elseif i == 12
                push!(row_data, "2103007000")  # SCC
            elseif i == 14
                push!(row_data, "100.0")  # ANN_VALUE
            elseif i == 18
                push!(row_data, "100.0")  # STKHGT
            elseif i == 19
                push!(row_data, "10.0")  # STKDIAM
            elseif i == 20
                push!(row_data, "500.0")  # STKTEMP
            elseif i == 21
                push!(row_data, "1000.0")  # STKFLOW
            elseif i == 22
                push!(row_data, "50.0")  # STKVEL
            elseif i in 53:64
                push!(row_data, "10.0")  # Monthly values
            else
                push!(row_data, "test")
            end
        end
        open(tmpfile, "w") do io
            println(io, join(row_data, ","))
        end
        result = read_ff10(tmpfile, :point)
        @test result isa FF10PointDataFrame
        rm(tmpfile)
    end

    @testset "read_ff10 invalid format" begin
        tmpfile = tempname() * ".csv"
        open(tmpfile, "w") do io
            println(io, "a,b,c")
        end
        @test_throws ArgumentError read_ff10(tmpfile, :invalid)
        rm(tmpfile)
    end

    @testset "read_ff10 with comments" begin
        # Create a temp CSV with 45 columns and comment lines
        tmpfile = tempname() * ".csv"
        row_data = []
        for i in 1:45
            if i == 2
                push!(row_data, "36001")  # FIPS
            elseif i == 6
                push!(row_data, "2103007000")  # SCC
            elseif i == 9
                push!(row_data, "100.0")  # ANN_VALUE
            elseif i in 21:32
                push!(row_data, "10.0")  # Monthly values
            else
                push!(row_data, "test")
            end
        end
        open(tmpfile, "w") do io
            println(io, "# This is a comment line")
            println(io, "# Another comment")
            println(io, join(row_data, ","))
            println(io, "# Comment in the middle")
            println(io, join(row_data, ","))
        end
        result = read_ff10(tmpfile, :nonpoint)
        @test result isa FF10NonPointDataFrame
        @test size(result.df, 1) == 2  # Only 2 data rows, comments ignored
        @test result.df[1, :FIPS] == "36001"
        @test result.df[2, :FIPS] == "36001"
        rm(tmpfile)
    end

    @testset "aggregate_emissions" begin
        df1 = DataFrame(
            POLID = ["NOX", "VOC"],
            COUNTRY = ["USA", "USA"],
            FIPS = ["36001", "36001"],
            SCC = ["2103007000", "2103007000"],
            ANN_VALUE = [100.0, 50.0]
        )
        df2 = DataFrame(
            POLID = ["NOX", "SO2"],
            COUNTRY = ["USA", "USA"],
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            ANN_VALUE = [200.0, 75.0]
        )
        result = aggregate_emissions([df1, df2])
        @test hasproperty(result, :LONGITUDE)
        @test hasproperty(result, :LATITUDE)
        # NOX for FIPS 36001 should be summed: 100 + 200 = 300
        nox_rows = filter(row -> row.POLID == "NOX" && row.FIPS == "36001", result)
        @test nrow(nox_rows) == 1
        @test nox_rows[1, :ANN_VALUE] == 300.0
        # Result should be sorted descending by ANN_VALUE
        @test result[1, :ANN_VALUE] >= result[end, :ANN_VALUE]
    end

    @testset "aggregate_emissions preserves existing coordinates" begin
        df_with_coords = DataFrame(
            POLID = ["NOX"],
            COUNTRY = ["USA"],
            FIPS = ["36001"],
            SCC = ["2103007000"],
            LONGITUDE = [-73.5],
            LATITUDE = [40.5],
            ANN_VALUE = [100.0]
        )
        result = aggregate_emissions([df_with_coords])
        @test result[1, :LONGITUDE] == -73.5
        @test result[1, :LATITUDE] == 40.5
    end

    @testset "filter_known_pollutants" begin
        df = DataFrame(
            POLID = ["NOX", "UNKNOWN_POLL", "VOC", "FAKE_POLL", "PM25-PRI"],
            ANN_VALUE = [100.0, 50.0, 75.0, 25.0, 60.0]
        )
        result = filter_known_pollutants(df)
        @test nrow(result) == 3
        @test "NOX" in result.POLID
        @test "VOC" in result.POLID
        @test "PM25-PRI" in result.POLID
        @test !("UNKNOWN_POLL" in result.POLID)
        @test !("FAKE_POLL" in result.POLID)
    end

    @testset "map_pollutant_names!" begin
        df = DataFrame(
            POLID = ["EXH__VOC", "PM25-PRI", "NOX", "SO2", "EVP__VOC"]
        )
        map_pollutant_names!(df)
        @test df.POLID[1] == "VOC"
        @test df.POLID[2] == "PM25"
        @test df.POLID[3] == "NOX"
        @test df.POLID[4] == "SO2"
        @test df.POLID[5] == "VOC"
    end

    @testset "map_pollutant_names! passthrough unknown" begin
        df = DataFrame(POLID = ["UNKNOWN_POLL"])
        map_pollutant_names!(df)
        @test df.POLID[1] == "UNKNOWN_POLL"
    end

    @testset "assign_surrogates direct match" begin
        emissions = DataFrame(
            POLID = ["NOX", "VOC"],
            COUNTRY = ["USA", "USA"],
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            ANN_VALUE = [100.0, 50.0]
        )
        gridref = DataFrame(
            COUNTRY = ["USA", "USA"],
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            Surrogate = [100, 200]
        )
        result = assign_surrogates(emissions, gridref)
        @test hasproperty(result, :Surrogate)
        @test result[1, :Surrogate] == 100
        @test result[2, :Surrogate] == 200
    end

    @testset "assign_surrogates fallback match" begin
        emissions = DataFrame(
            POLID = ["NOX"],
            COUNTRY = ["USA"],
            FIPS = ["36999"],  # Not in gridref
            SCC = ["2103007000"],
            ANN_VALUE = [100.0]
        )
        gridref = DataFrame(
            COUNTRY = ["USA", "USA"],
            FIPS = ["36001", "00000"],
            SCC = ["2103007000", "2103007000"],
            Surrogate = [100, 300]
        )
        result = assign_surrogates(emissions, gridref)
        @test result[1, :Surrogate] == 300  # Should get fallback surrogate
        @test result[1, :FIPS] == "36999"   # Original FIPS should be preserved
    end

    @testset "assign_surrogates no match" begin
        emissions = DataFrame(
            POLID = ["NOX"],
            COUNTRY = ["USA"],
            FIPS = ["36999"],
            SCC = ["9999999999"],
            ANN_VALUE = [100.0]
        )
        gridref = DataFrame(
            COUNTRY = ["USA"],
            FIPS = ["36001"],
            SCC = ["2103007000"],
            Surrogate = [100]
        )
        result = assign_surrogates(emissions, gridref)
        @test result[1, :Surrogate] == 0  # Unmatched records default to 0
    end

    @testset "build_data_weight_map" begin
        emissions = DataFrame(
            POLID = ["NOX", "VOC"],
            COUNTRY = ["USA", "USA"],
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            ANN_VALUE = [100.0, 50.0],
            Surrogate = [100, 200]
        )
        srgSpecs = [
            SurrogateSpec(
                "USA", "Pop", 100, "/data/pop.shp", "POP",
                "/weight/pop_weight.shp", "", String[], String[], Float64[], "", String[], Float64[]
            ),
            SurrogateSpec(
                "USA", "Area", 200, "/data/area.shp", "AREA",
                "/weight/area_weight.shp", "", String[], String[], Float64[], "", String[], Float64[]
            ),
            SurrogateSpec(
                "USA", "Unused", 999, "/data/unused.shp", "UNUSED",
                "/weight/unused_weight.shp", "", String[], String[], Float64[], "", String[], Float64[]
            ),
        ]
        result = build_data_weight_map(emissions, srgSpecs)
        @test length(result) == 2
        @test haskey(result, ("/data/pop.shp", "/weight/pop_weight.shp"))
        @test haskey(result, ("/data/area.shp", "/weight/area_weight.shp"))
        @test !haskey(result, ("/data/unused.shp", "/weight/unused_weight.shp"))
        @test "USA+100" in result[("/data/pop.shp", "/weight/pop_weight.shp")]
        @test "USA+200" in result[("/data/area.shp", "/weight/area_weight.shp")]
    end

    @testset "build_data_weight_map with missing surrogates" begin
        emissions = DataFrame(
            POLID = ["NOX"],
            COUNTRY = ["USA"],
            FIPS = ["36001"],
            SCC = ["2103007000"],
            ANN_VALUE = [100.0],
            Surrogate = [missing]
        )
        srgSpecs = SurrogateSpec[]
        result = build_data_weight_map(emissions, srgSpecs)
        @test length(result) == 0
    end

    @testset "build_data_weight_map without Surrogate column" begin
        emissions = DataFrame(
            POLID = ["NOX"],
            COUNTRY = ["USA"],
            FIPS = ["36001"],
            SCC = ["2103007000"],
            ANN_VALUE = [100.0]
        )
        srgSpecs = SurrogateSpec[]
        result = build_data_weight_map(emissions, srgSpecs)
        @test length(result) == 0
    end

    @testset "find_surrogate_by_code with region" begin
        srgs = [
            SurrogateSpec(
                "USA", "Pop", 100, "", "", "", "",
                String[], String[], Float64[], "", String[], Float64[]
            ),
            SurrogateSpec(
                "Canada", "Pop", 100, "", "", "", "",
                String[], String[], Float64[], "", String[], Float64[]
            ),
            SurrogateSpec(
                "USA", "Area", 200, "", "", "", "",
                String[], String[], Float64[], "", String[], Float64[]
            ),
        ]
        @test find_surrogate_by_code(srgs, "USA", 100).Region == "USA"
        @test find_surrogate_by_code(srgs, "Canada", 100).Region == "Canada"
        @test find_surrogate_by_code(srgs, "USA", 200).Name == "Area"
        @test find_surrogate_by_code(srgs, "Mexico", 100) === nothing
        @test find_surrogate_by_code(srgs, "USA", 999) === nothing
    end
end
