using LibGEOS, SparseArrays

@testset "Output tests" begin
    @testset "format_float" begin
        @test Emissions.format_float(0.0) == "0.0"
        @test Emissions.format_float(1e-25) == "0.0"
        @test Emissions.format_float(1.23456) == "1.234560"
        @test Emissions.format_float(1e-5) == "1.000000e-05"
        @test Emissions.format_float(1e7) == "1.000000e+07"
    end

    @testset "find_surrogate_by_code" begin
        srgs = [
            SurrogateSpec("US", "Pop", 100, "", "", "", "",
                String[], String[], Float64[], "", String[], Float64[]),
            SurrogateSpec("US", "Area", 200, "", "", "", "",
                String[], String[], Float64[], "", String[], Float64[]),
        ]
        @test find_surrogate_by_code(srgs, 100).Name == "Pop"
        @test find_surrogate_by_code(srgs, 200).Name == "Area"
        @test find_surrogate_by_code(srgs, 999) === nothing
    end

    @testset "get_data_weight_shapefiles" begin
        srg = SurrogateSpec("US", "Pop", 100, "/data.shp", "POP",
            "/weight.shp", "", String[], String[], Float64[], "", String[], Float64[])
        data_path, weight_path = get_data_weight_shapefiles(srg)
        @test data_path == "/data.shp"
        @test weight_path == "/weight.shp"
    end

    @testset "writeEmis" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        data = sparse([1, 2], [1, 2], [1.5, 2.5], 2, 2)

        tmpfile = tempname() * ".csv"
        writeEmis(tmpfile, data, grid; pollutant="NOX", units="kg/s")

        content = read(tmpfile, String)
        @test occursin("NOX", content)
        @test occursin("kg/s", content)
        @test occursin("1,1,", content)
        @test occursin("2,2,", content)
        rm(tmpfile)
    end
end
