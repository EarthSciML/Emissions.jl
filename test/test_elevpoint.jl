using DataFrames

@testset "Elevpoint tests" begin
    @testset "analytical_plume_rise basic" begin
        # Tall hot stack should have positive plume rise
        pr = analytical_plume_rise(100.0, 5.0, 450.0, 15.0)
        @test pr > 0.0

        # Zero velocity, minimal temperature difference
        pr = analytical_plume_rise(10.0, 1.0, 293.15, 0.0)
        @test pr ≈ 0.0

        # Cold stack with high velocity -> momentum dominated
        pr = analytical_plume_rise(10.0, 2.0, 300.0, 20.0)
        @test pr > 0.0
    end

    @testset "analytical_plume_rise increases with temperature" begin
        pr_low = analytical_plume_rise(50.0, 3.0, 350.0, 10.0)
        pr_high = analytical_plume_rise(50.0, 3.0, 500.0, 10.0)
        @test pr_high > pr_low
    end

    @testset "analytical_plume_rise increases with velocity" begin
        pr_low = analytical_plume_rise(50.0, 3.0, 400.0, 5.0)
        pr_high = analytical_plume_rise(50.0, 3.0, 400.0, 20.0)
        @test pr_high > pr_low
    end

    @testset "analytical_plume_rise increases with diameter" begin
        pr_small = analytical_plume_rise(50.0, 1.0, 400.0, 10.0)
        pr_large = analytical_plume_rise(50.0, 5.0, 400.0, 10.0)
        @test pr_large > pr_small
    end

    @testset "classify_point_sources surface sources" begin
        # Small stack below all thresholds
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1.0e-6],
            STKHGT = [5.0],     # 5m < 15.24m threshold
            STKDIAM = [0.5],
            STKTEMP = [300.0],  # 300K < 366.48K threshold
            STKVEL = [1.0],     # 1 m/s < 3.048 m/s threshold
        )
        result = classify_point_sources(df)
        @test result[1, :source_class] == "surface"
        @test hasproperty(result, :analytical_plume_rise)
    end

    @testset "classify_point_sources elevated by height" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1.0e-3],
            STKHGT = [20.0],    # > 15.24m
            STKDIAM = [1.0],
            STKTEMP = [300.0],
            STKVEL = [1.0],
        )
        result = classify_point_sources(df)
        @test result[1, :source_class] in ["elevated", "ping"]
    end

    @testset "classify_point_sources elevated by temperature" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1.0e-3],
            STKHGT = [5.0],
            STKDIAM = [1.0],
            STKTEMP = [400.0],  # > 366.48K
            STKVEL = [1.0],
        )
        result = classify_point_sources(df)
        @test result[1, :source_class] in ["elevated", "ping"]
    end

    @testset "classify_point_sources ping classification" begin
        # Tall stack with high emissions -> PinG
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1.0],
            STKHGT = [50.0],    # > 30.48m ping threshold
            STKDIAM = [3.0],
            STKTEMP = [450.0],
            STKVEL = [15.0],
        )
        result = classify_point_sources(df)
        @test result[1, :source_class] == "ping"
    end

    @testset "classify_point_sources mixed sources" begin
        df = DataFrame(
            FIPS = ["36001", "36001", "36001"],
            SCC = ["2103007000", "2103007000", "2103007000"],
            POLID = ["NOX", "NOX", "NOX"],
            ANN_VALUE = [1.0e-6, 1.0e-3, 1.0],
            STKHGT = [3.0, 20.0, 50.0],
            STKDIAM = [0.3, 1.0, 3.0],
            STKTEMP = [295.0, 350.0, 500.0],
            STKVEL = [0.5, 5.0, 20.0],
        )
        result = classify_point_sources(df)
        @test nrow(result) == 3
        @test result[1, :source_class] == "surface"
        # The other two should be at least elevated
        @test result[2, :source_class] in ["elevated", "ping"]
        @test result[3, :source_class] in ["elevated", "ping"]
    end

    @testset "classify_point_sources with custom criteria" begin
        criteria = ElevationCriteria(
            100.0,  # min_stack_height
            50.0,   # min_exit_velocity
            500.0,  # min_exit_temperature
            100.0,  # min_flow_rate
            100.0,  # min_plume_rise
            200.0,  # ping_stack_height
            0.0,    # ping_emissions_threshold
        )
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1.0],
            STKHGT = [20.0],   # < 100m custom threshold
            STKDIAM = [1.0],
            STKTEMP = [350.0], # < 500K custom threshold
            STKVEL = [5.0],    # < 50 m/s custom threshold
        )
        result = classify_point_sources(df; criteria = criteria)
        @test result[1, :source_class] == "surface"
    end

    @testset "group_stacks same facility" begin
        df = DataFrame(
            FIPS = ["36001", "36001", "36001"],
            STKHGT = [50.0, 55.0, 100.0],   # First two within 10m bin
            STKTEMP = [400.0, 420.0, 400.0], # First two within 50K bin
            STKDIAM = [3.0, 3.0, 5.0],
            STKVEL = [10.0, 10.0, 20.0],
            LONGITUDE = [-73.5, -73.5, -73.5],
            LATITUDE = [40.5, 40.5, 40.5],
        )
        result = group_stacks(df)
        @test hasproperty(result, :stack_group)
        # First two should be same group, third different
        @test result[1, :stack_group] == result[2, :stack_group]
        @test result[1, :stack_group] != result[3, :stack_group]
    end

    @testset "group_stacks different facilities" begin
        df = DataFrame(
            FIPS = ["36001", "36005"],
            STKHGT = [50.0, 50.0],
            STKTEMP = [400.0, 400.0],
            STKDIAM = [3.0, 3.0],
            STKVEL = [10.0, 10.0],
            LONGITUDE = [-73.5, -74.0],
            LATITUDE = [40.5, 41.0],
        )
        result = group_stacks(df)
        @test result[1, :stack_group] != result[2, :stack_group]
    end

    @testset "group_stacks custom bins" begin
        df = DataFrame(
            FIPS = ["36001", "36001"],
            STKHGT = [50.0, 70.0],     # 20m apart
            STKTEMP = [400.0, 410.0],
            STKDIAM = [3.0, 3.0],
            STKVEL = [10.0, 10.0],
            LONGITUDE = [-73.5, -73.5],
            LATITUDE = [40.5, 40.5],
        )
        # Default bin (10m) -> different groups
        result1 = group_stacks(df; height_bin = 10.0)
        @test result1[1, :stack_group] != result1[2, :stack_group]

        # Wider bin (25m) -> same group
        result2 = group_stacks(df; height_bin = 25.0)
        @test result2[1, :stack_group] == result2[2, :stack_group]
    end

    @testset "group_stacks without coordinates" begin
        # Without LONGITUDE/LATITUDE, groups by FIPS only
        df = DataFrame(
            FIPS = ["36001", "36001", "36005"],
            STKHGT = [50.0, 55.0, 50.0],
            STKTEMP = [400.0, 420.0, 400.0],
            STKDIAM = [3.0, 3.0, 3.0],
            STKVEL = [10.0, 10.0, 10.0],
        )
        result = group_stacks(df)
        @test result[1, :stack_group] == result[2, :stack_group]
        @test result[1, :stack_group] != result[3, :stack_group]
    end
end
