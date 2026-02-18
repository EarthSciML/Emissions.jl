using DataFrames, Dates, SparseArrays

@testset "Merge tests" begin
    @testset "merge_emissions basic" begin
        hourly = DataFrame(
            FIPS = ["36001", "36001"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)],
            emission_rate = [100.0, 200.0]
        )

        locIndex = Dict{String, IndexInfo}(
            "36001" => IndexInfo([1], [1], [1.0], true, true)
        )
        grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        result = merge_emissions(hourly, locIndex, grid)
        @test nrow(result) == 2  # 2 hours
        @test all(result.grid_row .== 1)
        @test all(result.grid_col .== 1)
        sort!(result, :hour)
        @test result[1, :emission_rate] ≈ 100.0
        @test result[2, :emission_rate] ≈ 200.0
    end

    @testset "merge_emissions distributes across cells" begin
        hourly = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            hour = [DateTime(2019, 7, 1, 0)],
            emission_rate = [100.0]
        )

        locIndex = Dict{String, IndexInfo}(
            "36001" => IndexInfo([1, 1], [1, 2], [0.6, 0.4], true, true)
        )
        grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        result = merge_emissions(hourly, locIndex, grid)
        @test nrow(result) == 2
        sort!(result, :grid_col)
        @test result[1, :emission_rate] ≈ 60.0
        @test result[2, :emission_rate] ≈ 40.0
        @test sum(result.emission_rate) ≈ 100.0
    end

    @testset "merge_emissions filters pollutants" begin
        hourly = DataFrame(
            FIPS = ["36001", "36001"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "UNKNOWN_POLL"],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            emission_rate = [100.0, 50.0]
        )

        locIndex = Dict{String, IndexInfo}(
            "36001" => IndexInfo([1], [1], [1.0], true, true)
        )
        grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        result = merge_emissions(hourly, locIndex, grid)
        @test nrow(result) == 1
        @test result[1, :pollutant] == "NOX"
    end

    @testset "merge_emissions skips missing locations" begin
        hourly = DataFrame(
            FIPS = ["99999"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            hour = [DateTime(2019, 7, 1, 0)],
            emission_rate = [100.0]
        )

        locIndex = Dict{String, IndexInfo}()
        grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        result = merge_emissions(hourly, locIndex, grid)
        @test nrow(result) == 0
    end

    @testset "merge_emissions accumulates same cell+hour+pollutant" begin
        hourly = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            emission_rate = [100.0, 50.0]
        )

        # Both FIPS map to the same grid cell
        locIndex = Dict{String, IndexInfo}(
            "36001" => IndexInfo([1], [1], [1.0], true, true),
            "36005" => IndexInfo([1], [1], [1.0], true, true)
        )
        grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        result = merge_emissions(hourly, locIndex, grid)
        @test nrow(result) == 1
        @test result[1, :emission_rate] ≈ 150.0  # 100 + 50
    end

    @testset "merge_emissions with point sources" begin
        hourly = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            hour = [DateTime(2019, 7, 1, 0)],
            emission_rate = [100.0],
            LONGITUDE = [-73.5],
            LATITUDE = [40.5]
        )

        key = location_key("36001", -73.5, 40.5)
        locIndex = Dict{String, IndexInfo}(
            key => IndexInfo([1], [1], [1.0], true, true)
        )
        grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        result = merge_emissions(hourly, locIndex, grid)
        @test nrow(result) == 1
        @test result[1, :emission_rate] ≈ 100.0
    end

    @testset "merge_categories basic" begin
        df1 = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 2],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "NOX"],
            emission_rate = [60.0, 40.0]
        )
        df2 = DataFrame(
            grid_row = [1],
            grid_col = [1],
            hour = [DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX"],
            emission_rate = [25.0]
        )

        result = merge_categories(df1, df2)
        @test nrow(result) == 2
        sort!(result, :grid_col)
        @test result[1, :emission_rate] ≈ 85.0  # 60 + 25
        @test result[2, :emission_rate] ≈ 40.0
    end

    @testset "merge_categories empty input" begin
        result = merge_categories(DataFrame[])
        @test nrow(result) == 0
        @test hasproperty(result, :emission_rate)
    end

    @testset "merge_categories single DataFrame" begin
        df = DataFrame(
            grid_row = [1],
            grid_col = [1],
            hour = [DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX"],
            emission_rate = [100.0]
        )
        result = merge_categories(df)
        @test nrow(result) == 1
        @test result[1, :emission_rate] ≈ 100.0
    end

    @testset "merge_categories preserves multiple pollutants" begin
        df1 = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "VOC"],
            emission_rate = [100.0, 50.0]
        )
        df2 = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "SO2"],
            emission_rate = [25.0, 75.0]
        )

        result = merge_categories(df1, df2)
        @test nrow(result) == 3  # NOX (summed), VOC, SO2
        nox = filter(r -> r.pollutant == "NOX", result)
        @test nrow(nox) == 1
        @test nox[1, :emission_rate] ≈ 125.0
    end

    @testset "merge_categories_tracked" begin
        df1 = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 2],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "NOX"],
            emission_rate = [60.0, 40.0]
        )
        df2 = DataFrame(
            grid_row = [1],
            grid_col = [1],
            hour = [DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX"],
            emission_rate = [25.0]
        )

        result = merge_categories_tracked(
            [
                "area" => df1,
                "point" => df2,
            ]
        )
        @test nrow(result) == 3
        @test hasproperty(result, :source_category)
        area = filter(r -> r.source_category == "area", result)
        @test nrow(area) == 2
        point = filter(r -> r.source_category == "point", result)
        @test nrow(point) == 1
        @test point[1, :emission_rate] ≈ 25.0
    end

    @testset "merge_categories_tracked empty" begin
        result = merge_categories_tracked(Pair{String, DataFrame}[])
        @test nrow(result) == 0
        @test hasproperty(result, :source_category)
    end

    @testset "merge_2d_3d" begin
        surface = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 2],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "NOX"],
            emission_rate = [100.0, 50.0]
        )
        elevated = DataFrame(
            grid_row = [1],
            grid_col = [1],
            hour = [DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX"],
            layer = [3],
            emission_rate = [25.0]
        )

        result = merge_2d_3d(surface, elevated)
        @test nrow(result) == 3
        @test hasproperty(result, :layer)
        # Surface should be layer 1
        surf_rows = filter(r -> r.layer == 1, result)
        @test nrow(surf_rows) == 2
        # Elevated should be layer 3
        elev_rows = filter(r -> r.layer == 3, result)
        @test nrow(elev_rows) == 1
        @test elev_rows[1, :emission_rate] ≈ 25.0
    end

    @testset "to_model_ready" begin
        merged = DataFrame(
            grid_row = [1, 1, 2],
            grid_col = [1, 2, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)],
            pollutant = ["NOX", "NOX", "VOC"],
            emission_rate = [100.0, 50.0, 75.0]
        )
        grid = NewGridRegular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        hours = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)]

        result = to_model_ready(merged, grid, hours)
        @test haskey(result, "NOX")
        @test haskey(result, "VOC")
        @test size(result["NOX"]) == (3, 3, 1, 2)
        @test result["NOX"][1, 1, 1, 1] ≈ 100.0
        @test result["NOX"][1, 2, 1, 1] ≈ 50.0
        @test result["VOC"][2, 1, 1, 2] ≈ 75.0
    end

    @testset "to_model_ready with layers" begin
        merged = DataFrame(
            grid_row = [1, 1],
            grid_col = [1, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0)],
            pollutant = ["NOX", "NOX"],
            layer = [1, 3],
            emission_rate = [100.0, 50.0]
        )
        grid = NewGridRegular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        hours = [DateTime(2019, 7, 1, 0)]

        result = to_model_ready(merged, grid, hours; n_layers = 5)
        @test size(result["NOX"]) == (2, 2, 5, 1)
        @test result["NOX"][1, 1, 1, 1] ≈ 100.0
        @test result["NOX"][1, 1, 3, 1] ≈ 50.0
        @test result["NOX"][1, 1, 2, 1] ≈ 0.0  # Unused layer
    end
end
