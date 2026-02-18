using Test
using Emissions
using DataFrames
using Dates

@testset "Pipeline Integration Tests" begin
    @testset "speciation → temporal → merge pipeline" begin
        # --- Setup: synthetic emissions with both point and area sources ---
        emissions = DataFrame(
            FIPS = ["36001", "36001", "36005", "36005"],
            SCC = ["2103007000", "2103007000", "2103007000", "2103007000"],
            POLID = ["NOX", "VOC", "NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0, 200.0, 80.0],
            COUNTRY = ["USA", "USA", "USA", "USA"],
            LONGITUDE = [-73.5, -73.5, missing, missing],
            LATITUDE = [40.5, 40.5, missing, missing],
        )

        # Speciation profiles: NOX → NO (90%) + NO2 (10%), VOC passed through
        gspro = DataFrame(
            profile_code = ["P001", "P001"],
            pollutant_id = ["NOX", "NOX"],
            species_id = ["NO", "NO2"],
            split_factor = [0.9, 0.1],
            divisor = [1.0, 1.0],
            mass_fraction = [0.9, 0.1],
        )
        gsref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            pollutant_id = ["NOX"],
            profile_code = ["P001"],
        )

        # Temporal profiles: uniform
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
            profile_id = [1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), fill(1.0 / 24.0, 24)],
        )
        xref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1],
        )

        # --- Step 1: Speciation ---
        speciated = speciate_emissions(emissions, gspro, gsref; basis = :mass)

        # Verify column naming: should have :POLID, not :species
        @test hasproperty(speciated, :POLID)
        @test !hasproperty(speciated, :species)

        # Verify location columns are preserved
        @test hasproperty(speciated, :LONGITUDE)
        @test hasproperty(speciated, :LATITUDE)
        @test hasproperty(speciated, :COUNTRY)

        # NOX sources should be split into NO and NO2
        # VOC sources pass through (no profile match)
        no_rows = filter(r -> r.POLID == "NO", speciated)
        no2_rows = filter(r -> r.POLID == "NO2", speciated)
        voc_rows = filter(r -> r.POLID == "VOC", speciated)
        @test nrow(no_rows) == 2   # one per FIPS
        @test nrow(no2_rows) == 2
        @test nrow(voc_rows) == 2  # passed through unchanged

        # Mass conservation for NOX: total speciated = total input
        total_nox_input = 100.0 + 200.0
        total_speciated_nox = sum(no_rows.ANN_VALUE) + sum(no2_rows.ANN_VALUE)
        @test total_speciated_nox ≈ total_nox_input

        # VOC mass unchanged
        @test sum(voc_rows.ANN_VALUE) ≈ 50.0 + 80.0

        # Verify SCC is normalized to 10 digits
        @test all(length.(speciated.SCC) .== 10)

        # Point source coordinates preserved for FIPS 36001
        point_rows = filter(r -> r.FIPS == "36001", speciated)
        @test all(r -> !ismissing(r.LONGITUDE), eachrow(point_rows))
        @test all(r -> r.LONGITUDE ≈ -73.5, eachrow(point_rows))
        @test all(r -> r.LATITUDE ≈ 40.5, eachrow(point_rows))

        # Area source coordinates remain missing for FIPS 36005
        area_rows = filter(r -> r.FIPS == "36005", speciated)
        @test all(r -> ismissing(r.LONGITUDE), eachrow(area_rows))

        # --- Step 2: Temporal allocation ---
        ep_start = DateTime(2019, 7, 1, 0)
        ep_end = DateTime(2019, 7, 1, 3)  # 3 hours

        hourly = temporal_allocate(speciated, profiles, xref, ep_start, ep_end)

        # Verify basic structure
        @test nrow(hourly) == nrow(speciated) * 3  # 3 hours per source
        @test hasproperty(hourly, :POLID)
        @test hasproperty(hourly, :hour)
        @test hasproperty(hourly, :emission_rate)

        # Verify location columns are preserved through temporal allocation
        @test hasproperty(hourly, :LONGITUDE)
        @test hasproperty(hourly, :LATITUDE)
        @test hasproperty(hourly, :COUNTRY)

        # Point source coordinates preserved
        hourly_point = filter(r -> r.FIPS == "36001", hourly)
        @test all(r -> !ismissing(r.LONGITUDE), eachrow(hourly_point))

        # Area source coordinates remain missing
        hourly_area = filter(r -> r.FIPS == "36005", hourly)
        @test all(r -> ismissing(r.LONGITUDE), eachrow(hourly_area))

        # With uniform profiles, hourly rate = annual rate (rate factor = 1.0)
        expected_rate_factor = 1.0
        for polid in unique(hourly.POLID)
            pol_speciated = filter(r -> r.POLID == polid, speciated)
            pol_hourly = filter(r -> r.POLID == polid, hourly)
            for fips in unique(pol_hourly.FIPS)
                src_ann = filter(r -> r.FIPS == fips, pol_speciated)
                src_hourly = filter(r -> r.FIPS == fips, pol_hourly)
                @test nrow(src_ann) == 1
                expected_rate = src_ann[1, :ANN_VALUE] * expected_rate_factor
                @test all(r -> isapprox(r, expected_rate, rtol = 1.0e-6), src_hourly.emission_rate)
            end
        end

        # --- Step 3: Spatial indexing ---
        # Use a simple lon/lat grid centered around the point source location
        grid = NewGridRegular("test", 3, 3, "EPSG:4326", 1.0, 1.0, -75.0, 39.0)

        locIndex = compute_grid_indices(hourly, grid)

        # Point source key should exist and be in grid
        point_key = location_key("36001", -73.5, 40.5)
        @test haskey(locIndex, point_key)
        @test locIndex[point_key].inGrid

        # Area source key should exist (though not in grid without shapefile)
        area_key = location_key("36005", missing, missing)
        @test haskey(locIndex, area_key)

        # --- Step 4: Merge emissions ---
        # Only point sources will have valid grid indices, so filter to those
        # Use species names as pollutant_groups for merge
        pol_list = unique(string.(hourly.POLID))
        merged = merge_emissions(hourly, locIndex, grid; species_list = pol_list)

        # Merged should have gridded hourly emissions for point sources at minimum
        @test nrow(merged) > 0
        @test hasproperty(merged, :grid_row)
        @test hasproperty(merged, :grid_col)
        @test hasproperty(merged, :hour)
        @test hasproperty(merged, :pollutant)
        @test hasproperty(merged, :emission_rate)

        # Verify all expected pollutants appear in merged output
        merged_pols = unique(merged.pollutant)
        @test "NO" in merged_pols
        @test "NO2" in merged_pols
        @test "VOC" in merged_pols

        # Mass conservation for point sources through merge:
        # For each pollutant, total merged rate should equal the hourly rate * fraction
        # Point sources are allocated entirely to one cell (fraction = 1.0)
        for pol in ["NO", "NO2", "VOC"]
            pol_hourly_point = filter(r -> r.POLID == pol && r.FIPS == "36001", hourly)
            pol_merged = filter(r -> r.pollutant == pol, merged)
            # Total merged rate should equal the hourly point-source rate per hour
            for hr in unique(pol_merged.hour)
                hr_merged = filter(r -> r.hour == hr, pol_merged)
                hr_hourly = filter(r -> r.hour == hr, pol_hourly_point)
                @test sum(hr_merged.emission_rate) ≈ sum(hr_hourly.emission_rate) rtol = 1.0e-6
            end
        end
    end

    @testset "temporal_allocate preserves location columns" begin
        emissions = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            ANN_VALUE = [100.0, 200.0],
            LONGITUDE = [-73.5, missing],
            LATITUDE = [40.5, missing],
            COUNTRY = ["USA", "USA"],
            Surrogate = [missing, 100],
        )
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
            profile_id = [1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), fill(1.0 / 24.0, 24)],
        )
        xref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1],
        )

        result = temporal_allocate(
            emissions, profiles, xref,
            DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 2)
        )

        @test hasproperty(result, :LONGITUDE)
        @test hasproperty(result, :LATITUDE)
        @test hasproperty(result, :COUNTRY)
        @test hasproperty(result, :Surrogate)
        @test nrow(result) == 4  # 2 sources * 2 hours

        # Verify values are correct per source
        point_rows = filter(r -> r.FIPS == "36001", result)
        @test all(r -> r.LONGITUDE ≈ -73.5, eachrow(point_rows))
        @test all(r -> r.LATITUDE ≈ 40.5, eachrow(point_rows))
        @test all(r -> r.COUNTRY == "USA", eachrow(point_rows))

        area_rows = filter(r -> r.FIPS == "36005", result)
        @test all(r -> ismissing(r.LONGITUDE), eachrow(area_rows))
        @test all(r -> r.Surrogate == 100, eachrow(area_rows))
    end

    @testset "temporal_allocate empty input with location columns" begin
        emissions = DataFrame(
            FIPS = String[],
            SCC = String[],
            POLID = String[],
            ANN_VALUE = Float64[],
            LONGITUDE = Float64[],
            LATITUDE = Float64[],
        )
        profiles = DataFrame(
            profile_type = String[],
            profile_id = Int[],
            factors = Vector{Float64}[],
        )
        xref = DataFrame(
            FIPS = String[],
            SCC = String[],
            monthly_id = Int[],
            weekly_id = Int[],
            diurnal_id = Int[],
        )

        result = temporal_allocate(
            emissions, profiles, xref,
            DateTime(2019, 1, 1), DateTime(2019, 1, 2)
        )
        @test nrow(result) == 0
        @test hasproperty(result, :LONGITUDE)
        @test hasproperty(result, :LATITUDE)
    end

    @testset "assign_surrogates coalesces missing values" begin
        emissions = DataFrame(
            POLID = ["NOX", "VOC"],
            COUNTRY = ["USA", "USA"],
            FIPS = ["36999", "36001"],
            SCC = ["9999999999", "2103007000"],
            ANN_VALUE = [100.0, 50.0],
        )
        gridref = DataFrame(
            COUNTRY = ["USA"],
            FIPS = ["36001"],
            SCC = ["2103007000"],
            Surrogate = [100],
        )

        result = assign_surrogates(emissions, gridref)
        # Matched row should have surrogate 100
        matched = filter(r -> r.FIPS == "36001", result)
        @test matched[1, :Surrogate] == 100

        # Unmatched row should have surrogate 0 (not missing)
        unmatched = filter(r -> r.FIPS == "36999", result)
        @test unmatched[1, :Surrogate] == 0
        @test !ismissing(unmatched[1, :Surrogate])
    end
end
