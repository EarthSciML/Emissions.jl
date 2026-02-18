using NCDatasets
using Dates

@testset "IOAPI Output" begin

    @testset "Write and read roundtrip" begin
        grid = NewGridRegular(
            "TEST_GRID", 4, 3, "+proj=longlat +datum=WGS84",
            1.0, 1.0, -100.0, 30.0
        )
        hours = [DateTime(2019, 7, 1, h) for h in 0:2]

        model_data = Dict{String, Array{Float64, 4}}()
        model_data["NO2"] = zeros(3, 4, 1, 3)
        model_data["NO2"][1, 1, 1, 1] = 1.5
        model_data["NO2"][2, 3, 1, 2] = 2.5
        model_data["SO2"] = ones(3, 4, 1, 3) * 0.1

        outfile = tempname() * ".nc"
        try
            write_ioapi(
                outfile, model_data, grid, hours;
                description = "Test file", history = "test run"
            )

            ds = NCDatasets.Dataset(outfile, "r")
            try
                # Check dimensions
                @test ds.dim["TSTEP"] == 3
                @test ds.dim["LAY"] == 1
                @test ds.dim["ROW"] == 3
                @test ds.dim["COL"] == 4
                @test ds.dim["VAR"] == 2
                @test ds.dim["DATE-TIME"] == 2

                # Check global attributes
                @test ds.attrib["GDTYP"] == Int32(1)  # LATGRD3
                @test ds.attrib["NCOLS"] == Int32(4)
                @test ds.attrib["NROWS"] == Int32(3)
                @test ds.attrib["NLAYS"] == Int32(1)
                @test ds.attrib["NVARS"] == Int32(2)
                @test ds.attrib["XORIG"] ≈ -100.0
                @test ds.attrib["YORIG"] ≈ 30.0
                @test ds.attrib["XCELL"] ≈ 1.0
                @test ds.attrib["YCELL"] ≈ 1.0
                @test ds.attrib["FTYPE"] == Int32(1)
                @test ds.attrib["FILEDESC"] == "Test file"
                @test ds.attrib["HISTORY"] == "test run"

                # Check TFLAG
                tflag = ds["TFLAG"][:, :, :]
                @test tflag[1, 1, 1] == Int32(2019182)  # July 1, 2019 = day 182
                @test tflag[2, 1, 1] == Int32(0)         # hour 0
                @test tflag[2, 1, 2] == Int32(10000)     # hour 1
                @test tflag[2, 1, 3] == Int32(20000)     # hour 2

                # Check variable data
                no2_data = ds["NO2"][:, :, :, :]
                @test no2_data[1, 1, 1, 1] ≈ Float32(1.5)
                @test no2_data[3, 2, 1, 2] ≈ Float32(2.5)
                @test no2_data[2, 2, 1, 1] ≈ Float32(0.0)

                so2_data = ds["SO2"][:, :, :, :]
                @test all(so2_data .≈ Float32(0.1))
            finally
                close(ds)
            end
        finally
            rm(outfile; force = true)
        end
    end

    @testset "SDATE and STIME encoding" begin
        grid = NewGridRegular(
            "TEST", 2, 2, "+proj=longlat +datum=WGS84",
            1.0, 1.0, 0.0, 0.0
        )
        hours = [DateTime(2020, 1, 15, 14, 30, 0)]

        model_data = Dict("CO" => ones(2, 2, 1, 1))

        outfile = tempname() * ".nc"
        try
            write_ioapi(outfile, model_data, grid, hours)
            ds = NCDatasets.Dataset(outfile, "r")
            try
                # Jan 15 = day 15
                @test ds.attrib["SDATE"] == Int32(2020015)
                @test ds.attrib["STIME"] == Int32(143000)
            finally
                close(ds)
            end
        finally
            rm(outfile; force = true)
        end
    end

    @testset "Multi-layer output" begin
        grid = NewGridRegular(
            "TEST", 3, 3, "+proj=longlat +datum=WGS84",
            1.0, 1.0, 0.0, 0.0
        )
        hours = [DateTime(2019, 7, 1)]

        model_data = Dict{String, Array{Float64, 4}}()
        model_data["PM25"] = zeros(3, 3, 3, 1)
        model_data["PM25"][1, 1, 1, 1] = 10.0  # layer 1
        model_data["PM25"][1, 1, 2, 1] = 5.0   # layer 2
        model_data["PM25"][1, 1, 3, 1] = 1.0   # layer 3

        outfile = tempname() * ".nc"
        try
            write_ioapi(outfile, model_data, grid, hours; n_layers = 3)
            ds = NCDatasets.Dataset(outfile, "r")
            try
                @test ds.dim["LAY"] == 3
                @test ds.attrib["NLAYS"] == Int32(3)
                pm_data = ds["PM25"][:, :, :, :]
                @test pm_data[1, 1, 1, 1] ≈ Float32(10.0)
                @test pm_data[1, 1, 2, 1] ≈ Float32(5.0)
                @test pm_data[1, 1, 3, 1] ≈ Float32(1.0)
            finally
                close(ds)
            end
        finally
            rm(outfile; force = true)
        end
    end

    @testset "LCC projection attributes" begin
        grid = NewGridRegular(
            "LCC_GRID", 10, 10,
            "+proj=lcc +lat_1=33.0 +lat_2=45.0 +lon_0=-97.0 +lat_0=40.0 +datum=WGS84",
            36000.0, 36000.0, -2736000.0, -2088000.0
        )
        hours = [DateTime(2019, 7, 1)]
        model_data = Dict("NOX" => ones(10, 10, 1, 1))

        outfile = tempname() * ".nc"
        try
            write_ioapi(outfile, model_data, grid, hours)
            ds = NCDatasets.Dataset(outfile, "r")
            try
                @test ds.attrib["GDTYP"] == Int32(2)  # LAMGRD3
                @test ds.attrib["P_ALP"] ≈ 33.0
                @test ds.attrib["P_BET"] ≈ 45.0
                @test ds.attrib["P_GAM"] ≈ -97.0
                @test ds.attrib["XCENT"] ≈ -97.0
                @test ds.attrib["YCENT"] ≈ 40.0
                @test ds.attrib["XCELL"] ≈ 36000.0
                @test ds.attrib["YCELL"] ≈ 36000.0
            finally
                close(ds)
            end
        finally
            rm(outfile; force = true)
        end
    end

    @testset "VAR-LIST attribute" begin
        grid = NewGridRegular(
            "TEST", 2, 2, "+proj=longlat +datum=WGS84",
            1.0, 1.0, 0.0, 0.0
        )
        hours = [DateTime(2019, 7, 1)]
        model_data = Dict(
            "CO" => ones(2, 2, 1, 1),
            "NO2" => ones(2, 2, 1, 1),
            "SO2" => ones(2, 2, 1, 1)
        )

        outfile = tempname() * ".nc"
        try
            write_ioapi(outfile, model_data, grid, hours)
            ds = NCDatasets.Dataset(outfile, "r")
            try
                varlist = ds.attrib["VAR-LIST"]
                # Variables are sorted alphabetically, each 16 chars
                @test length(varlist) == 3 * 16
                @test strip(varlist[1:16]) == "CO"
                @test strip(varlist[17:32]) == "NO2"
                @test strip(varlist[33:48]) == "SO2"
            finally
                close(ds)
            end
        finally
            rm(outfile; force = true)
        end
    end

    @testset "Convenience DataFrame method" begin
        grid = NewGridRegular(
            "TEST", 3, 3, "+proj=longlat +datum=WGS84",
            1.0, 1.0, 0.0, 0.0
        )
        hours = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)]

        merged = DataFrame(
            grid_row = [1, 1, 2],
            grid_col = [1, 2, 1],
            hour = [DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)],
            pollutant = ["CO", "CO", "CO"],
            emission_rate = [1.0, 2.0, 3.0]
        )

        outfile = tempname() * ".nc"
        try
            write_ioapi(outfile, merged, grid, hours)
            ds = NCDatasets.Dataset(outfile, "r")
            try
                @test ds.dim["TSTEP"] == 2
                @test ds.dim["ROW"] == 3
                @test ds.dim["COL"] == 3
                co_data = ds["CO"][:, :, :, :]
                @test co_data[1, 1, 1, 1] ≈ Float32(1.0)
                @test co_data[2, 1, 1, 1] ≈ Float32(2.0)
                @test co_data[1, 2, 1, 2] ≈ Float32(3.0)
            finally
                close(ds)
            end
        finally
            rm(outfile; force = true)
        end
    end

    @testset "Empty model_data error" begin
        grid = NewGridRegular(
            "TEST", 2, 2, "+proj=longlat +datum=WGS84",
            1.0, 1.0, 0.0, 0.0
        )
        hours = [DateTime(2019, 7, 1)]
        @test_throws ErrorException write_ioapi(tempname(), Dict{String, Array{Float64, 4}}(), grid, hours)
    end

    @testset "_parse_proj4_to_ioapi" begin
        # Lat-lon
        attrs = Emissions._parse_proj4_to_ioapi("+proj=longlat +datum=WGS84")
        @test attrs.GDTYP == Int32(1)

        # LCC
        attrs = Emissions._parse_proj4_to_ioapi("+proj=lcc +lat_1=33.0 +lat_2=45.0 +lon_0=-97.0 +lat_0=40.0")
        @test attrs.GDTYP == Int32(2)
        @test attrs.P_ALP ≈ 33.0
        @test attrs.P_BET ≈ 45.0
        @test attrs.P_GAM ≈ -97.0
        @test attrs.YCENT ≈ 40.0

        # UTM
        attrs = Emissions._parse_proj4_to_ioapi("+proj=utm +zone=17")
        @test attrs.GDTYP == Int32(5)
        @test attrs.P_ALP ≈ 17.0

        # Stereographic
        attrs = Emissions._parse_proj4_to_ioapi("+proj=stere +lat_0=90.0 +lon_0=-80.0")
        @test attrs.GDTYP == Int32(6)
        @test attrs.P_ALP ≈ 90.0
        @test attrs.P_GAM ≈ -80.0

        # Mercator
        attrs = Emissions._parse_proj4_to_ioapi("+proj=merc +lon_0=0.0")
        @test attrs.GDTYP == Int32(7)

        # Transverse Mercator
        attrs = Emissions._parse_proj4_to_ioapi("+proj=tmerc +lon_0=-75.0 +lat_0=40.0")
        @test attrs.GDTYP == Int32(8)
        @test attrs.P_GAM ≈ -75.0
        @test attrs.YCENT ≈ 40.0
    end

    @testset "_datetime_to_julian" begin
        jdate, jtime = Emissions._datetime_to_julian(DateTime(2019, 7, 1, 0, 0, 0))
        @test jdate == Int32(2019182)
        @test jtime == Int32(0)

        jdate, jtime = Emissions._datetime_to_julian(DateTime(2020, 1, 1, 12, 30, 45))
        @test jdate == Int32(2020001)
        @test jtime == Int32(123045)

        # Leap year: Dec 31, 2020 = day 366
        jdate, jtime = Emissions._datetime_to_julian(DateTime(2020, 12, 31, 23, 59, 59))
        @test jdate == Int32(2020366)
        @test jtime == Int32(235959)
    end

    @testset "Roundtrip coordtype <-> proj4" begin
        # LCC roundtrip
        proj4 = Emissions._coordtype_to_proj4(2, 33.0, 45.0, -97.0, -97.0, 40.0)
        attrs = Emissions._parse_proj4_to_ioapi(proj4)
        @test attrs.GDTYP == Int32(2)
        @test attrs.P_ALP ≈ 33.0
        @test attrs.P_BET ≈ 45.0
        @test attrs.P_GAM ≈ -97.0
        @test attrs.YCENT ≈ 40.0

        # Lat-lon roundtrip
        proj4 = Emissions._coordtype_to_proj4(1, 0.0, 0.0, 0.0, 0.0, 0.0)
        attrs = Emissions._parse_proj4_to_ioapi(proj4)
        @test attrs.GDTYP == Int32(1)
    end
end
