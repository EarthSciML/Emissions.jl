using DataFrames
import GeoInterface as GI
import GeometryOps as GO

@testset "Types tests" begin
    @testset "SurrogateSpec construction" begin
        srg = SurrogateSpec(
            "US", "Population", 100, "/path/data.shp", "POP",
            "/path/weight.shp", "Population surrogate",
            String["Backup1"], String["POP2020"], Float64[1.0],
            "", String[], Float64[]
        )
        @test srg.Region == "US"
        @test srg.Name == "Population"
        @test srg.Code == 100
        @test srg.DataShapefile == "/path/data.shp"
        @test srg.DataAttribute == "POP"
    end

    @testset "GridDef construction" begin
        grid = GridDef(
            "test", 1, 1, "EPSG:4326",
            [[(0.0, 0.0), (1.0, 1.0)]]
        )
        @test grid.Name == "test"
        @test grid.Nx == 1
        @test grid.Ny == 1
    end

    @testset "Config construction" begin
        cfg = Config(
            ["gridref1.txt"],
            "srgspec.csv",
            "/shapefiles",
            "+proj=longlat",
            "+proj=lcc",
            "grid.txt",
            "TestGrid",
            "counties.shp",
            "/output"
        )
        @test cfg.GridName == "TestGrid"
        @test length(cfg.f_gridRef) == 1
    end

    @testset "IndexInfo construction" begin
        idx = IndexInfo([1, 2], [3, 4], [0.6, 0.4], true, true)
        @test idx.rows == [1, 2]
        @test idx.cols == [3, 4]
        @test sum(idx.fracs) ≈ 1.0
        @test idx.inGrid == true
        @test idx.coveredByGrid == true
    end

    @testset "SpatialProcessor construction" begin
        srg = SurrogateSpec(
            "US", "Pop", 100, "", "", "", "",
            String[], String[], Float64[], "", String[], Float64[]
        )
        grid = GridDef(
            "test", 1, 1, "EPSG:4326",
            [[(0.0, 0.0), (1.0, 1.0)]]
        )
        gridRef = DataFrame(SCC = String[], SrgCode = Int[])
        sp = SpatialProcessor([srg], grid, gridRef, "+proj=longlat", false, 100, 10)
        @test sp.MatchFullSCC == false
        @test sp.MemCacheSize == 100
    end
end
