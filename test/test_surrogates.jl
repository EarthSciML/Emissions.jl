using SparseArrays

@testset "Surrogates tests" begin
    @testset "generate_countySurrogate" begin
        data = Dict(
            "01001" => sparse([1, 2], [1, 1], [0.5, 0.5], 2, 2),
            "01002" => sparse([1], [2], [1.0], 2, 2),
        )
        weight = sparse([1, 2], [1, 2], [2.0, 3.0], 2, 2)
        grid = sparse([1, 1, 2, 2], [1, 2, 1, 2], [1.0, 1.0, 1.0, 1.0], 2, 2)

        result = generate_countySurrogate(data, weight, grid)

        @test haskey(result, "01001")
        @test haskey(result, "01002")

        # Check that fractions sum to 1 for county 01001
        @test sum(result["01001"]) ≈ 1.0
    end

    @testset "generate_countySurrogate zero weight fallback" begin
        data = Dict(
            "01001" => sparse([1], [1], [1.0], 2, 2),
        )
        weight = spzeros(2, 2)
        grid = sparse([1, 1, 2, 2], [1, 2, 1, 2], [1.0, 1.0, 1.0, 1.0], 2, 2)

        result = generate_countySurrogate(data, weight, grid)
        @test sum(result["01001"]) ≈ 1.0
    end

    @testset "update_locIndex with matching FIPS" begin
        locIndex = Dict{String, IndexInfo}()
        surrogates = Dict(
            "01001" => sparse([1, 2], [1, 2], [0.6, 0.4], 3, 3),
        )

        idx = update_locIndex(locIndex, "01001", surrogates)
        @test idx.inGrid == true
        @test length(idx.rows) == 2
        @test haskey(locIndex, "01001")
    end

    @testset "update_locIndex without matching FIPS" begin
        locIndex = Dict{String, IndexInfo}()
        surrogates = Dict{String, SparseMatrixCSC{Float64, Int}}()

        idx = update_locIndex(locIndex, "99999", surrogates)
        @test idx.inGrid == false
        @test isempty(idx.rows)
    end

    @testset "find_column_name" begin
        using DataFrames
        df = DataFrame(Population = [1, 2], Area = [3.0, 4.0])
        @test Emissions.find_column_name(df, "population") == "Population"
        @test Emissions.find_column_name(df, "POPULATION") == "Population"
        @test_throws ErrorException Emissions.find_column_name(df, "nonexistent")
    end

    @testset "read_crs_epsg" begin
        # Test with known CRS string
        test_crs = "GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]"

        # Create temporary .prj file
        temp_prj = tempname() * ".prj"
        open(temp_prj, "w") do io
            write(io, test_crs)
        end

        result = Emissions.read_crs_epsg(temp_prj)
        @test result isa Int
        @test result == 4326  # WGS84 EPSG code

        # Test with non-existent file
        @test_throws ErrorException Emissions.read_crs_epsg("nonexistent.prj")

        # Cleanup
        rm(temp_prj)
    end

    @testset "generate_data_sparse_matrices" begin
        # Create a small test grid
        test_grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Test error handling with non-existent shapefile
        @test_throws Exception generate_data_sparse_matrices(
            "nonexistent.shp", "ATTR", test_grid, "+proj=longlat +datum=WGS84"
        )

        # Test with empty attribute name
        @test_throws Exception generate_data_sparse_matrices(
            "nonexistent.shp", "", test_grid, "+proj=longlat +datum=WGS84"
        )
    end

    @testset "generate_weight_sparse_matrices" begin
        # Create a small test grid
        test_grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Test error handling with non-existent shapefile
        @test_throws Exception generate_weight_sparse_matrices(
            "nonexistent.shp", ["WEIGHT"], [1.0], test_grid, "+proj=longlat +datum=WGS84"
        )

        # Test with empty weight columns
        @test_throws Exception generate_weight_sparse_matrices(
            "nonexistent.shp", String[], Float64[], test_grid, "+proj=longlat +datum=WGS84"
        )
    end

    @testset "generate_grid_sparse_matrices" begin
        # Create a small test grid
        test_grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Test basic grid matrix generation
        result = generate_grid_sparse_matrices(test_grid)
        @test result isa SparseArrays.AbstractSparseArray
        # The actual matrix dimensions depend on the implementation
        @test size(result, 1) >= 1  # At least one row
        @test size(result, 2) >= 1  # At least one column

        # Test with single cell grid
        single_grid = NewGridIrregular("single", 1, 1, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        result2 = generate_grid_sparse_matrices(single_grid)
        @test result2 isa SparseArrays.AbstractSparseArray
        @test size(result2, 1) == 1
    end
end
