using DataFrames
import GeoInterface as GI
import GeometryOps as GO

@testset "Regridding tests" begin
    @testset "grid_polygons" begin
        grid = NewGridIrregular("test", 3, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        polys = grid_polygons(grid)
        @test length(polys) == 6  # 3x2 grid
        # First cell should be at (0,0)-(1,1)
        @test GO.area(polys[1]) ≈ 1.0
    end

    @testset "build_regridder" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Source geometry covering the entire grid
        src = [GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])]
        regridder = build_regridder(src, grid)

        # The regridder should have 4 destination cells and 1 source polygon
        @test size(regridder.intersections) == (4, 1)

        # Each cell should have non-zero intersection
        for i in 1:4
            @test regridder.intersections[i, 1] > 0.0
        end
    end

    @testset "build_regridder partial overlap" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Source geometry covering only the first cell
        src = [GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])]
        regridder = build_regridder(src, grid)

        # Only cell 1 should have intersection
        @test regridder.intersections[1, 1] > 0.0
    end

    @testset "build_regridder multiple sources" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        src = [
            GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]]),
            GI.Polygon([[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]]),
        ]
        regridder = build_regridder(src, grid)

        @test size(regridder.intersections) == (4, 2)
    end

    @testset "regrid mass conservation" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Source covering full grid
        src = [GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])]
        regridder = build_regridder(src, grid)

        # Total intersection area should equal source area
        total_intersection = sum(regridder.intersections[:, 1])
        source_area = GO.area(src[1])
        @test total_intersection ≈ source_area rtol = 1.0e-6
    end
end
