using LibGEOS, SparseArrays

@testset "Spatial tests" begin
    @testset "NewPolygon" begin
        coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        poly = NewPolygon(coords)
        @test LibGEOS.area(poly) ≈ 1.0
    end

    @testset "NewGridIrregular" begin
        grid = NewGridIrregular("test", 3, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        @test grid.Nx == 3
        @test grid.Ny == 2
        @test length(grid.Cells) == 6
        @test LibGEOS.area(grid.Cells[1]) ≈ 1.0
    end

    @testset "GetIndex point" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Point inside grid
        idx = GetIndex(0.5, 0.5, grid)
        @test idx.inGrid == true
        @test length(idx.rows) == 1

        # Point outside grid
        idx2 = GetIndex(5.0, 5.0, grid)
        @test idx2.inGrid == false
    end

    @testset "GetIndex polygon" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        # Polygon covering part of the grid
        poly = NewPolygon([(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5)])
        idx = GetIndex(poly, grid)
        @test idx.inGrid == true
        @test length(idx.rows) >= 1
    end

    @testset "recordToGrid" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        idx = IndexInfo([1, 2], [1, 2], [0.6, 0.4], true, true)
        result = recordToGrid(100.0, idx, 2, 2)
        @test result[1, 1] ≈ 60.0
        @test result[2, 2] ≈ 40.0
    end

    @testset "GridFactors" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        areas = GridFactors(grid)
        @test all(areas .≈ 1.0)
    end

    @testset "uniqueCoordinates" begin
        lons = [1.0, 2.0, 1.0, 3.0]
        lats = [4.0, 5.0, 4.0, 6.0]
        indices = uniqueCoordinates(lons, lats)
        @test length(indices) == 3
        @test indices == [1, 2, 4]
    end

    @testset "uniqueLoc" begin
        lons = [1.0, 2.0, 1.0, 3.0]
        lats = [4.0, 5.0, 4.0, 6.0]
        result = uniqueLoc(lons, lats)
        @test length(result) == 3
        @test result[(1.0, 4.0)] == 1
        @test result[(2.0, 5.0)] == 2
        @test result[(3.0, 6.0)] == 4
    end
end
