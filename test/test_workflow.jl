using DataFrames, SparseArrays, Unitful

@testset "Workflow tests" begin
    @testset "location_key" begin
        # Point source with coordinates
        @test location_key("36001", -73.5, 40.5) == "36001_-73.500000_40.500000"

        # Area source without coordinates
        @test location_key("36001", missing, missing) == "36001"

        # Missing one coordinate treated as area source
        @test location_key("36001", -73.5, missing) == "36001"
        @test location_key("36001", missing, 40.5) == "36001"

        # Integer coordinates
        @test location_key("00001", 1, 2) == "00001_1.000000_2.000000"
    end

    @testset "compute_grid_indices point source inside grid" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(
            FIPS = ["00001"],
            LONGITUDE = [0.5],
            LATITUDE = [0.5],
        )
        locIndex = compute_grid_indices(emissions, grid)
        key = location_key("00001", 0.5, 0.5)
        @test haskey(locIndex, key)
        @test locIndex[key].inGrid == true
        @test length(locIndex[key].rows) == 1
        @test locIndex[key].fracs[1] == 1.0
    end

    @testset "compute_grid_indices point source outside grid" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(
            FIPS = ["00002"],
            LONGITUDE = [10.0],
            LATITUDE = [10.0],
        )
        locIndex = compute_grid_indices(emissions, grid)
        key = location_key("00002", 10.0, 10.0)
        @test haskey(locIndex, key)
        @test locIndex[key].inGrid == false
    end

    @testset "compute_grid_indices area source without shapefile" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(FIPS = ["36001"])
        locIndex = compute_grid_indices(emissions, grid)
        @test haskey(locIndex, "36001")
        @test locIndex["36001"].inGrid == false  # No counties shapefile
    end

    @testset "compute_grid_indices deduplication" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(
            FIPS = ["00001", "00001", "00001"],
            LONGITUDE = [0.5, 0.5, 1.5],
            LATITUDE = [0.5, 0.5, 1.5],
        )
        locIndex = compute_grid_indices(emissions, grid)
        # Two unique coordinate pairs: (0.5, 0.5) and (1.5, 1.5)
        @test length(locIndex) == 2
    end

    @testset "compute_grid_indices mixed point and area" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(
            FIPS = ["00001", "36001"],
            LONGITUDE = [0.5, missing],
            LATITUDE = [0.5, missing],
        )
        locIndex = compute_grid_indices(emissions, grid)
        # Point source key and area source key
        @test length(locIndex) == 2
        point_key = location_key("00001", 0.5, 0.5)
        @test locIndex[point_key].inGrid == true
        @test locIndex["36001"].inGrid == false  # No shapefile
    end

    @testset "refine_indices_with_surrogates refines area sources" begin
        locIndex = Dict{String, IndexInfo}(
            "01001" => IndexInfo([1], [1], [1.0], true, true),
        )

        county_surrogates = Dict(
            "01001" => sparse([1, 1], [1, 2], [0.7, 0.3], 2, 2),
        )

        result = refine_indices_with_surrogates(locIndex, county_surrogates)
        @test length(result["01001"].rows) == 2
        @test sum(result["01001"].fracs) ≈ 1.0
    end

    @testset "refine_indices_with_surrogates leaves point sources unchanged" begin
        locIndex = Dict{String, IndexInfo}(
            "01001" => IndexInfo([1], [1], [1.0], true, true),
            "01001_-73.500000_40.500000" => IndexInfo([2], [2], [1.0], true, true),
        )

        county_surrogates = Dict(
            "01001" => sparse([1, 1], [1, 2], [0.7, 0.3], 2, 2),
        )

        result = refine_indices_with_surrogates(locIndex, county_surrogates)

        # Area source should be refined
        @test length(result["01001"].rows) == 2

        # Point source should be unchanged
        @test result["01001_-73.500000_40.500000"].rows == [2]
        @test result["01001_-73.500000_40.500000"].fracs == [1.0]
    end

    @testset "refine_indices_with_surrogates no matching surrogate" begin
        locIndex = Dict{String, IndexInfo}(
            "99999" => IndexInfo([1], [1], [1.0], true, true),
        )
        county_surrogates = Dict{String, SparseMatrixCSC{Float64, Int}}()

        result = refine_indices_with_surrogates(locIndex, county_surrogates)
        # Should remain unchanged since no surrogate matches
        @test result["99999"].rows == [1]
        @test result["99999"].fracs == [1.0]
    end

    @testset "allocate_emissions_to_grid single record" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        locIndex = Dict{String, IndexInfo}(
            "00001_0.500000_0.500000" => IndexInfo([1], [1], [1.0], true, true),
        )

        emissions = DataFrame(
            FIPS = ["00001"],
            POLID = ["NOX"],
            SCC = ["2103007000"],
            ANN_VALUE = [1.0e-3],
            LONGITUDE = [0.5],
            LATITUDE = [0.5],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 1
        @test result[1, :cellIndex] == 1
        @test result[1, :NOX] ≈ 1.0e-3
        @test result[1, :VOC] ≈ 0.0
        @test result[1, :SCC] == "2103007000"
    end

    @testset "allocate_emissions_to_grid multiple pollutants same cell" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        locIndex = Dict{String, IndexInfo}(
            "00001_0.500000_0.500000" => IndexInfo([1], [1], [1.0], true, true),
        )

        emissions = DataFrame(
            FIPS = ["00001", "00001"],
            POLID = ["NOX", "VOC"],
            SCC = ["2103007000", "2103007000"],
            ANN_VALUE = [1.0e-3, 5.0e-4],
            LONGITUDE = [0.5, 0.5],
            LATITUDE = [0.5, 0.5],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 1  # Same cell and SCC
        @test result[1, :NOX] ≈ 1.0e-3
        @test result[1, :VOC] ≈ 5.0e-4
    end

    @testset "allocate_emissions_to_grid distributes across multiple cells" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        locIndex = Dict{String, IndexInfo}(
            "00001" => IndexInfo([1, 1], [1, 2], [0.6, 0.4], true, true),
        )

        emissions = DataFrame(
            FIPS = ["00001"],
            POLID = ["NOX"],
            SCC = ["2103007000"],
            ANN_VALUE = [100.0],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 2  # Two cells
        sort!(result, :cellIndex)
        @test result[1, :NOX] ≈ 60.0   # 100 * 0.6
        @test result[2, :NOX] ≈ 40.0   # 100 * 0.4
    end

    @testset "allocate_emissions_to_grid skips unknown pollutants" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        locIndex = Dict{String, IndexInfo}(
            "00001" => IndexInfo([1], [1], [1.0], true, true),
        )

        emissions = DataFrame(
            FIPS = ["00001"],
            POLID = ["UNKNOWN_POLL"],
            SCC = ["2103007000"],
            ANN_VALUE = [100.0],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 0
    end

    @testset "allocate_emissions_to_grid skips locations outside grid" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        locIndex = Dict{String, IndexInfo}(
            "00001_10.000000_10.000000" => IndexInfo(Int[], Int[], Float64[], false, false),
        )

        emissions = DataFrame(
            FIPS = ["00001"],
            POLID = ["NOX"],
            SCC = ["2103007000"],
            ANN_VALUE = [100.0],
            LONGITUDE = [10.0],
            LATITUDE = [10.0],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 0
    end

    @testset "allocate_emissions_to_grid with Unitful values" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        locIndex = Dict{String, IndexInfo}(
            "00001_0.500000_0.500000" => IndexInfo([1], [1], [1.0], true, true),
        )

        emissions = DataFrame(
            FIPS = ["00001"],
            POLID = ["NOX"],
            SCC = ["2103007000"],
            ANN_VALUE = [1.0e-3u"kg/s"],
            LONGITUDE = [0.5],
            LATITUDE = [0.5],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 1
        @test result[1, :NOX] ≈ 1.0e-3
    end

    @testset "allocate_emissions_to_grid empty DataFrame" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        locIndex = Dict{String, IndexInfo}()

        emissions = DataFrame(
            FIPS = String[],
            POLID = String[],
            SCC = String[],
            ANN_VALUE = Float64[],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 0
        @test hasproperty(result, :cellIndex)
        @test hasproperty(result, :NOX)
    end

    @testset "allocate_emissions_to_grid groups by SCC" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        locIndex = Dict{String, IndexInfo}(
            "00001_0.500000_0.500000" => IndexInfo([1], [1], [1.0], true, true),
        )

        emissions = DataFrame(
            FIPS = ["00001", "00001"],
            POLID = ["NOX", "NOX"],
            SCC = ["2103007000", "9999999999"],  # Different SCCs
            ANN_VALUE = [100.0, 50.0],
            LONGITUDE = [0.5, 0.5],
            LATITUDE = [0.5, 0.5],
        )

        result = allocate_emissions_to_grid(emissions, locIndex, grid)
        @test nrow(result) == 2  # Two different SCCs in same cell
        @test sum(result.NOX) ≈ 150.0
    end

    @testset "process_emissions_spatial with point sources" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(
            FIPS = ["00001", "00001", "00002"],
            POLID = ["NOX", "VOC", "NOX"],
            SCC = ["2103007000", "2103007000", "2103007000"],
            ANN_VALUE = [1.0e-3, 5.0e-4, 2.0e-3],
            LONGITUDE = [0.5, 0.5, 1.5],
            LATITUDE = [0.5, 0.5, 0.5],
        )

        result = process_emissions_spatial(emissions, grid)
        @test nrow(result) >= 1
        total_nox = sum(result.NOX)
        @test total_nox ≈ 3.0e-3  # 1.0e-3 + 2.0e-3
        total_voc = sum(result.VOC)
        @test total_voc ≈ 5.0e-4
    end

    @testset "process_emissions_spatial with surrogates" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(
            FIPS = ["01001"],
            POLID = ["NOX"],
            SCC = ["2103007000"],
            ANN_VALUE = [100.0],
        )

        county_surrogates = Dict(
            "01001" => sparse([1, 1], [1, 2], [0.7, 0.3], 2, 2),
        )

        result = process_emissions_spatial(
            emissions, grid;
            county_surrogates = county_surrogates
        )
        @test nrow(result) == 2
        sort!(result, :cellIndex)
        @test result[1, :NOX] ≈ 70.0   # 100 * 0.7
        @test result[2, :NOX] ≈ 30.0   # 100 * 0.3
    end

    @testset "process_emissions_spatial empty input" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        emissions = DataFrame(
            FIPS = String[],
            POLID = String[],
            SCC = String[],
            ANN_VALUE = Float64[],
        )

        result = process_emissions_spatial(emissions, grid)
        @test nrow(result) == 0
    end

    @testset "End-to-end workflow with pipeline functions" begin
        # Step 1: Create synthetic emissions (pre-aggregated, no FF10 wrapper needed)
        emissions_data = DataFrame(
            POLID = ["NOX", "VOC", "NOX", "SO2", "PM25-PRI", "UNKNOWN"],
            COUNTRY = ["USA", "USA", "USA", "USA", "USA", "USA"],
            FIPS = ["00001", "00001", "00002", "00002", "00001", "00001"],
            SCC = ["2103007000", "2103007000", "2103007000", "2103007000", "2103007000", "2103007000"],
            ANN_VALUE = [150.0, 75.0, 200.0, 125.0, 50.0, 999.0],
            LONGITUDE = [0.5, 0.5, 1.5, 1.5, 0.5, 0.5],
            LATITUDE = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
        )

        # Step 2: Filter to known pollutants
        filtered = filter_known_pollutants(emissions_data)
        @test nrow(filtered) == 5  # UNKNOWN filtered out

        # Step 3: Map pollutant names
        map_pollutant_names!(filtered)
        @test "PM25" in filtered.POLID  # PM25-PRI mapped to PM25

        # Step 4: Spatially allocate to grid
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        result = process_emissions_spatial(filtered, grid)

        # Verify results
        @test nrow(result) >= 1
        @test hasproperty(result, :cellIndex)
        @test hasproperty(result, :NOX)
        @test hasproperty(result, :VOC)
        @test hasproperty(result, :SO2)
        @test hasproperty(result, :PM25)
        @test hasproperty(result, :SCC)

        # Total emissions should be preserved
        @test sum(result.NOX) ≈ 350.0   # 150 + 200
        @test sum(result.VOC) ≈ 75.0
        @test sum(result.SO2) ≈ 125.0
        @test sum(result.PM25) ≈ 50.0
    end
end
