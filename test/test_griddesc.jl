@testset "GRIDDESC Reader" begin

    @testset "Lambert Conformal Conic grid" begin
        griddesc_content = """
        ' '
        'LCC_US'
        2 33.0 45.0 -97.0 -97.0 40.0
        ' '
        'US36_CRO'
        'LCC_US' -2736000.0 -2088000.0 36000.0 36000.0 148 112 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            grid = read_griddesc(path, "US36_CRO")
            @test grid.Name == "US36_CRO"
            @test grid.Nx == 148
            @test grid.Ny == 112
            @test grid.Dx == 36000.0
            @test grid.Dy == 36000.0
            @test grid.X0 == -2736000.0
            @test grid.Y0 == -2088000.0
            @test is_regular(grid)
            @test occursin("+proj=lcc", grid.SR)
            @test occursin("+lat_1=33.0", grid.SR)
            @test occursin("+lat_2=45.0", grid.SR)
            @test occursin("+lon_0=-97.0", grid.SR)
            @test occursin("+lat_0=40.0", grid.SR)
            @test length(grid.Extent) == 148 * 112
        finally
            rm(path; force = true)
        end
    end

    @testset "Lat-lon grid" begin
        griddesc_content = """
        ' '
        'LATLON'
        1 0.0 0.0 0.0 0.0 0.0
        ' '
        'GLOBAL_1DEG'
        'LATLON' -180.0 -90.0 1.0 1.0 360 180 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            grid = read_griddesc(path, "GLOBAL_1DEG")
            @test grid.Name == "GLOBAL_1DEG"
            @test grid.Nx == 360
            @test grid.Ny == 180
            @test grid.Dx == 1.0
            @test grid.Dy == 1.0
            @test grid.X0 == -180.0
            @test grid.Y0 == -90.0
            @test occursin("+proj=longlat", grid.SR)
        finally
            rm(path; force = true)
        end
    end

    @testset "UTM grid" begin
        griddesc_content = """
        ' '
        'UTM17'
        5 17.0 0.0 0.0 0.0 0.0
        ' '
        'UTM_GRID'
        'UTM17' 500000.0 4000000.0 1000.0 1000.0 50 50 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            grid = read_griddesc(path, "UTM_GRID")
            @test grid.Name == "UTM_GRID"
            @test grid.Nx == 50
            @test grid.Ny == 50
            @test occursin("+proj=utm", grid.SR)
            @test occursin("+zone=17", grid.SR)
        finally
            rm(path; force = true)
        end
    end

    @testset "Multi-grid file" begin
        griddesc_content = """
        ' '
        'LCC_US'
        2 33.0 45.0 -97.0 -97.0 40.0
        ' '
        'CONUS_36KM'
        'LCC_US' -2736000.0 -2088000.0 36000.0 36000.0 148 112 1
        'CONUS_12KM'
        'LCC_US' -2556000.0 -1728000.0 12000.0 12000.0 396 246 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            grid36 = read_griddesc(path, "CONUS_36KM")
            @test grid36.Nx == 148
            @test grid36.Ny == 112
            @test grid36.Dx == 36000.0

            grid12 = read_griddesc(path, "CONUS_12KM")
            @test grid12.Nx == 396
            @test grid12.Ny == 246
            @test grid12.Dx == 12000.0
        finally
            rm(path; force = true)
        end
    end

    @testset "Grid not found error" begin
        griddesc_content = """
        ' '
        'LCC_US'
        2 33.0 45.0 -97.0 -97.0 40.0
        ' '
        'EXISTING_GRID'
        'LCC_US' -2736000.0 -2088000.0 36000.0 36000.0 148 112 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            @test_throws ErrorException read_griddesc(path, "NONEXISTENT_GRID")
        finally
            rm(path; force = true)
        end
    end

    @testset "Missing coordinate system error" begin
        griddesc_content = """
        ' '
        'LATLON'
        1 0.0 0.0 0.0 0.0 0.0
        ' '
        'BAD_GRID'
        'MISSING_COORD' -180.0 -90.0 1.0 1.0 360 180 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            @test_throws ErrorException read_griddesc(path, "BAD_GRID")
        finally
            rm(path; force = true)
        end
    end

    @testset "_coordtype_to_proj4" begin
        @test occursin("+proj=longlat", Emissions._coordtype_to_proj4(1, 0.0, 0.0, 0.0, 0.0, 0.0))
        @test occursin("+proj=lcc", Emissions._coordtype_to_proj4(2, 33.0, 45.0, -97.0, -97.0, 40.0))
        @test occursin("+proj=utm", Emissions._coordtype_to_proj4(5, 17.0, 0.0, 0.0, 0.0, 0.0))
        @test occursin("+proj=stere", Emissions._coordtype_to_proj4(6, 90.0, 0.0, -80.0, 0.0, 0.0))
        @test occursin("+proj=merc", Emissions._coordtype_to_proj4(7, 0.0, 0.0, 0.0, 0.0, 0.0))
        @test occursin("+proj=tmerc", Emissions._coordtype_to_proj4(8, 0.0, 0.0, -75.0, 0.0, 40.0))
        @test_throws ErrorException Emissions._coordtype_to_proj4(99, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    @testset "GRIDDESC with comment lines" begin
        griddesc_content = """
        # This is a comment
        ' '
        # Coordinate systems
        'LATLON'
        1 0.0 0.0 0.0 0.0 0.0
        ' '
        # Grid definitions
        'SIMPLE'
        'LATLON' -10.0 -10.0 0.5 0.5 40 40 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            grid = read_griddesc(path, "SIMPLE")
            @test grid.Nx == 40
            @test grid.Ny == 40
            @test grid.Dx == 0.5
        finally
            rm(path; force = true)
        end
    end

    @testset "Cell geometry from GRIDDESC grid" begin
        griddesc_content = """
        ' '
        'LATLON'
        1 0.0 0.0 0.0 0.0 0.0
        ' '
        'SMALL'
        'LATLON' 0.0 0.0 1.0 1.0 3 3 1
        ' '
        """
        path = tempname()
        write(path, griddesc_content)
        try
            grid = read_griddesc(path, "SMALL")
            # Check first cell bounds
            xmin, xmax, ymin, ymax = cell_bounds(grid, 1)
            @test xmin ≈ 0.0
            @test xmax ≈ 1.0
            @test ymin ≈ 0.0
            @test ymax ≈ 1.0

            # Check last cell (3,3) = index 9
            xmin, xmax, ymin, ymax = cell_bounds(grid, 9)
            @test xmin ≈ 2.0
            @test xmax ≈ 3.0
            @test ymin ≈ 2.0
            @test ymax ≈ 3.0

            # Point lookup
            idx = GetIndex(0.5, 0.5, grid)
            @test idx.inGrid
            @test idx.rows == [1]
            @test idx.cols == [1]
        finally
            rm(path; force = true)
        end
    end
end
