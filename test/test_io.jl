using DataFrames

@testset "IO tests" begin
    @testset "strip_missing" begin
        df = DataFrame(a=[1, missing, 3], b=["x", "y", missing])
        Emissions.strip_missing(df)
        @test df[2, :a] == ""
        @test df[3, :b] == ""
        @test df[1, :a] == 1
    end

    @testset "getCountry" begin
        @test Emissions.getCountry("01001") == "US"
        @test Emissions.getCountry("90001") == "US"
        @test Emissions.getCountry("10001") == "Mexico"
        @test Emissions.getCountry("20001") == "Canada"
        @test Emissions.getCountry("30001") == "Unknown"
        @test Emissions.getCountry("") == "Unknown"
    end

    @testset "read_grid" begin
        # Create a temp file with grid data
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Header comment")
            println(io, "001001, -85.5, 32.5")
            println(io, "001002, -86.0, 33.0")
        end
        data = Emissions.read_grid(tmpfile)
        @test size(data, 1) == 2
        @test data.FIPS[1] == "01001"
        @test data.FIPS[2] == "01002"
        @test data.Longitude[1] ≈ -85.5
        @test data.Latitude[2] ≈ 33.0
        rm(tmpfile)
    end

    @testset "getShapefilePath" begin
        # Create a temp directory with a shapefile
        dir = mktempdir()
        touch(joinpath(dir, "test_file.shp"))
        result = Emissions.getShapefilePath(dir, "test_file", false)
        @test endswith(result, "test_file.shp")

        result2 = Emissions.getShapefilePath(dir, "nonexistent", false)
        @test result2 == ""
        rm(dir, recursive=true)
    end

    @testset "getShapefilePath nonexistent dir" begin
        dir = mktempdir()
        result = Emissions.getShapefilePath(dir, "myfile", true)
        @test result == ""
        rm(dir, recursive=true)
    end

    @testset "validateShapefile" begin
        @test Emissions.validateShapefile("") == false
        tmpfile = tempname() * ".shp"
        touch(tmpfile)
        @test Emissions.validateShapefile(tmpfile) == true
        rm(tmpfile)
        @test Emissions.validateShapefile("/nonexistent/file.shp") == false
    end

    @testset "readSrgSpecSMOKE" begin
        tmpfile = tempname() * ".csv"
        open(tmpfile, "w") do io
            println(io, "# Comment line")
            println(io, "US,Population,100,pop.shp,POP,weight.shp,Population surrogate,,,,,,")
        end
        dir = mktempdir()
        srgs = Emissions.readSrgSpecSMOKE(tmpfile, dir, false)
        @test length(srgs) == 1
        @test srgs[1].Region == "US"
        @test srgs[1].Name == "Population"
        @test srgs[1].Code == 100
        rm(tmpfile)
        rm(dir, recursive=true)
    end
end
