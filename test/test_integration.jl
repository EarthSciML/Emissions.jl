using Test
using Emissions
using DataFrames
using CSV
using Unitful
using SparseArrays

@testset "Integration Tests" begin
    @testset "End-to-end emissions processing workflow" begin
        # Create synthetic FF10 data similar to tutorial
        synthetic_data = DataFrame(
            COUNTRY = ["0", "0", "0", "0"],
            FIPS = ["36001", "36001", "36005", "36005"],
            TRIBAL_CODE = ["0", "0", "0", "0"],
            CENSUS_TRACT = ["0", "0", "0", "0"],
            SHAPE_ID = ["0", "0", "0", "0"],
            SCC = ["2103007000", "2103007000", "2103007000", "2103007000"],
            EMIS_TYPE = ["", "", "", ""],
            POLID = ["NOX", "VOC", "NOX", "VOC"],
            ANN_VALUE = [150.5, 75.2, 200.1, 125.8],
            ANN_PCT_RED = [0.0, 0.0, 0.0, 0.0],
            CONTROL_IDS = ["", "", "", ""],
            CONTROL_MEASURES = ["", "", "", ""],
            CURRENT_COST = [0.0, 0.0, 0.0, 0.0],
            CUMULATIVE_COST = [0.0, 0.0, 0.0, 0.0],
            PROJECTION_FACTOR = [1.0, 1.0, 1.0, 1.0],
            REG_CODES = ["", "", "", ""],
            CALC_METHOD = [1, 1, 1, 1],
            CALC_YEAR = [2019, 2019, 2019, 2019],
            DATE_UPDATED = ["", "", "", ""],
            DATA_SET_ID = ["", "", "", ""],
            JAN_VALUE = [0.0, 0.0, 0.0, 0.0],
            FEB_VALUE = [0.0, 0.0, 0.0, 0.0],
            MAR_VALUE = [0.0, 0.0, 0.0, 0.0],
            APR_VALUE = [0.0, 0.0, 0.0, 0.0],
            MAY_VALUE = [0.0, 0.0, 0.0, 0.0],
            JUN_VALUE = [0.0, 0.0, 0.0, 0.0],
            JUL_VALUE = [0.0, 0.0, 0.0, 0.0],
            AUG_VALUE = [0.0, 0.0, 0.0, 0.0],
            SEP_VALUE = [0.0, 0.0, 0.0, 0.0],
            OCT_VALUE = [0.0, 0.0, 0.0, 0.0],
            NOV_VALUE = [0.0, 0.0, 0.0, 0.0],
            DEC_VALUE = [0.0, 0.0, 0.0, 0.0],
            JAN_PCTRED = [0.0, 0.0, 0.0, 0.0],
            FEB_PCTRED = [0.0, 0.0, 0.0, 0.0],
            MAR_PCTRED = [0.0, 0.0, 0.0, 0.0],
            APR_PCTRED = [0.0, 0.0, 0.0, 0.0],
            MAY_PCTRED = [0.0, 0.0, 0.0, 0.0],
            JUN_PCTRED = [0.0, 0.0, 0.0, 0.0],
            JUL_PCTRED = [0.0, 0.0, 0.0, 0.0],
            AUG_PCTRED = [0.0, 0.0, 0.0, 0.0],
            SEP_PCTRED = [0.0, 0.0, 0.0, 0.0],
            OCT_PCTRED = [0.0, 0.0, 0.0, 0.0],
            NOV_PCTRED = [0.0, 0.0, 0.0, 0.0],
            DEC_PCTRED = [0.0, 0.0, 0.0, 0.0],
            COMMENT = ["Test", "Test", "Test", "Test"]
        )

        # Save original value before FF10NonPointDataFrame mutates the DataFrame in-place
        original_value = synthetic_data[1, :ANN_VALUE]  # tons/year (plain Float64)

        # Test FF10 data loading and unit conversion
        ff10_data = FF10NonPointDataFrame(synthetic_data)
        processed_emis = ff10_data.df

        @test size(processed_emis, 1) == 4
        @test size(processed_emis, 2) == 45

        # Verify unit conversion happened (tons/year to kg/s)
        converted_value = processed_emis[1, :ANN_VALUE]  # kg/s (with units)
        # Check that units are present and the value was converted
        @test unit(converted_value) == u"kg/s"
        @test ustrip(converted_value) < original_value  # Should be much smaller after conversion (comparing numeric values)

        # Test aggregation and filtering
        grouped_emis = combine(
            groupby(processed_emis, [:POLID, :COUNTRY, :FIPS, :SCC]),
            :ANN_VALUE => sum => :ANN_VALUE
        )
        @test size(grouped_emis, 1) == 4  # Should have 4 unique combinations

        # Test pollutant mapping
        known_polls = filter(row -> haskey(Pollutants, row.POLID), grouped_emis)
        @test size(known_polls, 1) == 4  # All pollutants should be known

        known_polls[!, :POLID] = [Pollutants[p] for p in known_polls[!, :POLID]]
        @test "NOX" in known_polls[!, :POLID]  # The mapping keeps "NOX" as "NOX"
        @test "VOC" in known_polls[!, :POLID]
    end

    @testset "Spatial processing components integration" begin
        # Test grid creation
        test_grid = NewGridIrregular("TestGrid", 2, 2, "+proj=lcc", 1000.0, 1000.0, 0.0, 0.0)

        @test test_grid.Name == "TestGrid"
        @test test_grid.Nx == 2
        @test test_grid.Ny == 2
        @test length(test_grid.Cells) == 4
        @test length(test_grid.Extent) == 4

        # Test GetIndex function with point
        test_point_idx = GetIndex(-100.0, 40.0, test_grid)
        @test test_point_idx isa IndexInfo

        # Test GridFactors (takes only GridDef)
        factors = GridFactors(test_grid)
        @test factors isa Matrix{Float64}
        @test size(factors) == (2, 2)  # Should match grid dimensions

        # Test coordinate utilities
        lon_vec = [-100.0, -100.0, -101.0, -101.0]
        lat_vec = [40.0, 40.0, 41.0, 41.0]
        unique_coords = uniqueCoordinates(lon_vec, lat_vec)
        @test unique_coords isa Vector{Int}
        @test length(unique_coords) <= 4  # Should deduplicate

        # Test location mapping
        lons = [-100.0, -100.0, -101.0, -101.0]
        lats = [40.0, 40.0, 41.0, 41.0]
        location_map = uniqueLoc(lons, lats)
        @test location_map isa Dict
    end

    @testset "Configuration and setup integration" begin
        # Create temporary files for configuration test
        temp_gridref = tempname() * ".csv"
        temp_srgspec = tempname() * ".csv"
        temp_dir = mktempdir()

        # Create test grid for this test set
        test_grid = NewGridIrregular("TestGrid", 2, 2, "+proj=lcc", 1000.0, 1000.0, 0.0, 0.0)

        # Write test grid reference file
        test_gridref = DataFrame(
            COUNTRY = ["USA", "USA"],
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            Surrogate = [100, 100]
        )
        CSV.write(temp_gridref, test_gridref)

        # Write test surrogate specification file (without headers for SMOKE format)
        test_srgspec_data = DataFrame(
            col1 = ["USA"],
            col2 = ["Test Population"],
            col3 = [100],
            col4 = ["pop.shp"],
            col5 = ["POP"],
            col6 = ["area.shp"],
            col7 = ["Test surrogate"],
            col8 = [""],
            col9 = ["AREA"],
            col10 = [1.0],
            col11 = [""],
            col12 = [""],
            col13 = [""]
        )
        CSV.write(temp_srgspec, test_srgspec_data; header=false)

        # Create test_srgspec for later use
        test_srgspec = DataFrame(
            Region = ["USA"],
            Name = ["Test Population"],
            Code = [100],
            DataShapefile = ["pop.shp"],
            DataAttribute = ["POP"],
            WeightShapefile = ["area.shp"],
            Details = ["Test surrogate"],
            BackupSurrogateNames = [""],
            WeightColumns = ["AREA"],
            WeightFactors = [1.0],
            FilterFunction = [""],
            MergeNames = [""],
            MergeMultipliers = [missing]
        )

        # Create test configuration
        config = Config(
            [temp_gridref],
            temp_srgspec,
            temp_dir,
            "+proj=longlat +datum=WGS84",
            "+proj=lcc +lat_1=33 +lat_2=45",
            "/nonexistent/grid.txt",  # Non-existent grid file
            "TestGrid",
            "/nonexistent/counties.shp",
            "/tmp/output/"
        )

        # Test NewSpatialProcessor creation
        # Convert DataFrame row to SurrogateSpec
        row = test_srgspec[1, :]
        test_surrogatespec = SurrogateSpec(
            row.Region, row.Name, row.Code, row.DataShapefile, row.DataAttribute,
            row.WeightShapefile, row.Details, [row.BackupSurrogateNames],
            [row.WeightColumns], [row.WeightFactors], row.FilterFunction,
            [row.MergeNames], ismissing(row.MergeMultipliers) ? Float64[] : [row.MergeMultipliers]
        )
        sp = NewSpatialProcessor([test_surrogatespec], test_grid, test_gridref, config.InputSR, false)
        @test sp isa SpatialProcessor
        @test length(sp.SrgSpecs) == 1
        @test sp.GridRef isa DataFrame

        # Test setupSpatialProcessor (should handle missing grid file gracefully)
        @test_logs (:warn, r"Grid file .* not found") begin
            sp2 = setupSpatialProcessor(config)
            @test sp2 isa SpatialProcessor
            @test sp2.Grids isa GridDef
            @test sp2.Grids.Nx == 1  # Should fall back to 1x1 grid
            @test sp2.Grids.Ny == 1
        end

        # Cleanup
        rm(temp_gridref)
        rm(temp_srgspec)
        rm(temp_dir, recursive=true)
    end

    @testset "Unit conversion integration" begin
        # Test that all conversion constants are consistent
        @test tonperyear isa typeof(1.0u"kg/s")
        @test tonpermonth isa typeof(1.0u"kg/s")
        @test foot isa typeof(1.0u"m")

        # Test kelvin function
        temps_f = [32.0, 68.0, 212.0]
        for temp_f in temps_f
            temp_k = kelvin(temp_f)
            @test temp_k isa typeof(1.0u"K")
            @test ustrip(temp_k) > 0  # All should be positive Kelvin
        end

        # Test freezing point
        freezing_k = kelvin(32.0)
        @test isapprox(ustrip(freezing_k), 273.15, rtol=1e-10)

        # Test that conversion factors make sense
        @test ustrip(tonperyear) > 0
        @test ustrip(tonpermonth) > ustrip(tonperyear)  # Monthly should be larger rate
        @test isapprox(ustrip(foot), 0.3048, rtol=1e-10)  # feet to meters
    end

    @testset "Pollutant mapping integration" begin
        # Test that pollutant dictionary maps FF10 codes to standard names
        @test haskey(Pollutants, "NOX")
        @test haskey(Pollutants, "VOC")
        @test haskey(Pollutants, "PM25-PRI")  # Use actual key from dictionary
        @test haskey(Pollutants, "SO2")
        @test haskey(Pollutants, "NH3")

        # Test mappings are reasonable (note: many map to themselves)
        @test Pollutants["NOX"] == "NOX"
        @test Pollutants["VOC"] == "VOC"
        @test Pollutants["SO2"] == "SO2"
        @test Pollutants["NH3"] == "NH3"

        # Test that unknown pollutants can be filtered out
        test_df = DataFrame(POLID = ["NOX", "UNKNOWN_POLLUTANT", "VOC"])
        known = filter(row -> haskey(Pollutants, row.POLID), test_df)
        @test nrow(known) == 2  # Should keep only NOX and VOC
    end
end