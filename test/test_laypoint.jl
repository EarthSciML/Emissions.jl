using DataFrames

@testset "Laypoint tests" begin
    @testset "LayerConfig construction" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0])
        @test config.n_layers == 5
        @test config.layer_heights[1] == 0.0
        @test config.layer_heights[end] == 5000.0
    end

    @testset "LayerConfig requires at least 2 heights" begin
        @test_throws ArgumentError LayerConfig([0.0])
    end

    @testset "compute_layer_fractions sum to 1" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0])

        # Plume within first two layers
        fracs = compute_layer_fractions(50.0, 300.0, config)
        @test length(fracs) == 5
        @test sum(fracs) ≈ 1.0

        # Plume at surface
        fracs = compute_layer_fractions(0.0, 50.0, config)
        @test sum(fracs) ≈ 1.0

        # Plume across all layers
        fracs = compute_layer_fractions(0.0, 4000.0, config)
        @test sum(fracs) ≈ 1.0
    end

    @testset "compute_layer_fractions single layer" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0])

        # Plume entirely within first layer
        fracs = compute_layer_fractions(10.0, 50.0, config)
        @test fracs[1] ≈ 1.0
        @test fracs[2] ≈ 0.0
        @test fracs[3] ≈ 0.0
    end

    @testset "compute_layer_fractions point source" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0])

        # Point source (zero plume depth)
        fracs = compute_layer_fractions(250.0, 250.0, config)
        @test fracs[2] ≈ 1.0  # 250m is in layer 2 (100-500m)
    end

    @testset "compute_layer_fractions above model top" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0])

        # Plume extending above model top
        fracs = compute_layer_fractions(800.0, 1500.0, config)
        @test sum(fracs) ≈ 1.0
        @test fracs[3] > 0.0  # Top layer gets the excess
    end

    @testset "compute_layer_fractions proportional distribution" begin
        config = LayerConfig([0.0, 100.0, 200.0, 300.0])

        # Plume spanning exactly layers 1 and 2 equally
        fracs = compute_layer_fractions(50.0, 150.0, config)
        @test fracs[1] ≈ 0.5  # 50m overlap in 100m layer
        @test fracs[2] ≈ 0.5  # 50m overlap in 100m layer
        @test fracs[3] ≈ 0.0
    end

    @testset "allocate_point_to_layers basic" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0])
        met = MetProfile(
            fill(280.0, 5),    # temperature
            fill(5.0, 5),      # wind speed
            fill(0.0, 5),      # stability class (unstable)
            fill(0.0, 5),      # stability param
        )

        # Hot tall stack -> should have plume rise
        fracs = allocate_point_to_layers(100.0, 5.0, 500.0, 20.0, met, config)
        @test length(fracs) == 5
        @test sum(fracs) ≈ 1.0
    end

    @testset "allocate_point_to_layers surface source" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0])
        met = MetProfile(
            fill(280.0, 3),
            fill(5.0, 3),
            fill(0.0, 3),
            fill(0.0, 3),
        )

        # Very small stack -> mostly in first layer
        fracs = allocate_point_to_layers(5.0, 0.5, 295.0, 1.0, met, config)
        @test sum(fracs) ≈ 1.0
        @test fracs[1] > 0.5  # Most emissions near surface
    end

    @testset "laypoint basic operation" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0])
        met = MetProfile(
            fill(280.0, 3),
            fill(5.0, 3),
            fill(0.0, 3),
            fill(0.0, 3),
        )
        met_profiles = Dict("default" => met)

        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0],
            STKHGT = [50.0],
            STKDIAM = [3.0],
            STKTEMP = [450.0],
            STKVEL = [15.0],
        )

        result = laypoint(df, met_profiles, config)
        @test hasproperty(result, :layer)
        @test hasproperty(result, :layer_fraction)
        @test nrow(result) >= 1
        @test sum(result.layer_fraction) ≈ 1.0
        # All original columns preserved
        @test hasproperty(result, :FIPS)
        @test hasproperty(result, :ANN_VALUE)
    end

    @testset "laypoint mass conservation" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0])
        met = MetProfile(
            fill(280.0, 5),
            fill(5.0, 5),
            fill(0.0, 5),
            fill(0.0, 5),
        )
        met_profiles = Dict("default" => met)

        df = DataFrame(
            FIPS = ["36001", "36001"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0],
            STKHGT = [80.0, 200.0],
            STKDIAM = [3.0, 5.0],
            STKTEMP = [400.0, 500.0],
            STKVEL = [10.0, 20.0],
        )

        result = laypoint(df, met_profiles, config)

        # Total fraction per source should be 1.0
        for polid in ["NOX", "VOC"]
            source_rows = filter(r -> r.POLID == polid, result)
            @test sum(source_rows.layer_fraction) ≈ 1.0
        end
    end

    @testset "laypoint multiple sources" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0])
        met = MetProfile(
            fill(280.0, 3),
            fill(5.0, 3),
            fill(0.0, 3),
            fill(0.0, 3),
        )
        met_profiles = Dict("default" => met)

        df = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            ANN_VALUE = [100.0, 200.0],
            STKHGT = [20.0, 150.0],
            STKDIAM = [1.0, 5.0],
            STKTEMP = [350.0, 500.0],
            STKVEL = [5.0, 20.0],
        )

        result = laypoint(df, met_profiles, config)
        @test nrow(result) >= 2  # At least one row per source

        # Check each source sums to 1.0
        for fips in ["36001", "36005"]
            source_rows = filter(r -> r.FIPS == fips, result)
            @test sum(source_rows.layer_fraction) ≈ 1.0
        end
    end

    @testset "laypoint without met data" begin
        config = LayerConfig([0.0, 100.0, 500.0, 1000.0])
        met_profiles = Dict{String, MetProfile}()  # No met data

        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0],
            STKHGT = [50.0],
            STKDIAM = [3.0],
            STKTEMP = [400.0],
            STKVEL = [10.0],
        )

        result = laypoint(df, met_profiles, config)
        @test nrow(result) == 1
        @test result[1, :layer] == 1  # Surface layer
        @test result[1, :layer_fraction] ≈ 1.0
    end
end
