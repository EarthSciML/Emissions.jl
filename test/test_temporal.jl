using DataFrames, Dates, Unitful

@testset "Temporal tests" begin
    @testset "read_temporal_profiles" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Comment line")
            println(io, "/PROFILE/")
            println(io, "MONTHLY, 1, 0.08, 0.08, 0.08, 0.09, 0.09, 0.09, 0.09, 0.09, 0.08, 0.08, 0.08, 0.08")
            println(io, "WEEKLY, 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0")
            println(io, "DIURNAL, 1, " * join(fill("0.041667", 24), ", "))
        end
        profiles = read_temporal_profiles(tmpfile)
        @test nrow(profiles) == 3
        @test profiles[1, :profile_type] == "MONTHLY"
        @test profiles[1, :profile_id] == 1
        @test length(profiles[1, :factors]) == 12
        @test sum(profiles[1, :factors]) ≈ 1.0 atol = 0.01
        @test profiles[2, :profile_type] == "WEEKLY"
        @test length(profiles[2, :factors]) == 7
        @test sum(profiles[2, :factors]) ≈ 7.0
        @test profiles[3, :profile_type] == "DIURNAL"
        @test length(profiles[3, :factors]) == 24
        rm(tmpfile)
    end

    @testset "read_temporal_profiles empty file" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Only comments")
            println(io, "")
        end
        profiles = read_temporal_profiles(tmpfile)
        @test nrow(profiles) == 0
        rm(tmpfile)
    end

    @testset "read_temporal_xref" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Comment line")
            println(io, "036001;2103007000;1;2;3")
            println(io, "000000;2103007000;4;5;6!national default")
        end
        xref = read_temporal_xref(tmpfile)
        @test nrow(xref) == 2
        @test xref[1, :FIPS] == "36001"
        @test xref[1, :SCC] == "2103007000"
        @test xref[1, :monthly_id] == 1
        @test xref[1, :weekly_id] == 2
        @test xref[1, :diurnal_id] == 3
        @test xref[2, :FIPS] == "00000"
        @test xref[2, :monthly_id] == 4
        rm(tmpfile)
    end

    @testset "read_temporal_xref with 5-digit FIPS" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "36001;2103007000;1;2;3")
        end
        xref = read_temporal_xref(tmpfile)
        @test nrow(xref) == 1
        @test xref[1, :FIPS] == "36001"
        rm(tmpfile)
    end

    @testset "_match_temporal_xref hierarchical matching" begin
        xref = DataFrame(
            FIPS = ["36001", "00000", "36005"],
            SCC = ["2103007000", "2103007000", "2103007000"],
            monthly_id = [10, 20, 30],
            weekly_id = [11, 21, 31],
            diurnal_id = [12, 22, 32]
        )

        # Exact FIPS+SCC match
        result = Emissions._match_temporal_xref(xref, "36001", "2103007000")
        @test result.monthly_id == 10

        # Falls back to national default (FIPS=00000)
        result = Emissions._match_temporal_xref(xref, "99999", "2103007000")
        @test result.monthly_id == 20

        # No match at all -> default (1,1,1)
        result = Emissions._match_temporal_xref(xref, "99999", "9999999999")
        @test result.monthly_id == 1
        @test result.weekly_id == 1
        @test result.diurnal_id == 1
    end

    @testset "_lookup_profile found" begin
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY"],
            profile_id = [1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7)]
        )
        result = Emissions._lookup_profile(profiles, "MONTHLY", 1)
        @test length(result) == 12
        @test sum(result) ≈ 1.0
    end

    @testset "_lookup_profile not found returns uniform default" begin
        profiles = DataFrame(
            profile_type = String[],
            profile_id = Int[],
            factors = Vector{Float64}[]
        )
        monthly = Emissions._lookup_profile(profiles, "MONTHLY", 999)
        @test length(monthly) == 12
        @test sum(monthly) ≈ 1.0

        weekly = Emissions._lookup_profile(profiles, "WEEKLY", 999)
        @test length(weekly) == 7
        @test sum(weekly) ≈ 7.0

        diurnal = Emissions._lookup_profile(profiles, "DIURNAL", 999)
        @test length(diurnal) == 24
        @test sum(diurnal) ≈ 1.0
    end

    @testset "temporal_allocate uniform profiles" begin
        # With uniform profiles, hourly rate should equal annual rate
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0]
        )
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
            profile_id = [1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), fill(1.0 / 24.0, 24)]
        )
        xref = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        ep_start = DateTime(2019, 7, 1, 0)
        ep_end = DateTime(2019, 7, 1, 3)  # 3 hours

        result = temporal_allocate(emissions, profiles, xref, ep_start, ep_end)

        @test nrow(result) == 3  # 3 hours
        @test all(result.POLID .== "NOX")
        @test all(result.FIPS .== "36001")
        # With uniform profiles: rate = 100 * (1/12) * (1/7) * (1/24 * 24) = 100/84
        # Actually: rate = ann * mf * (wf/7) * (df * 24)
        # With uniform: mf=1/12, wf=1.0, df=1/24
        # rate = 100 * (1/12) * (1/7) * (1/24 * 24) = 100 * (1/12) * (1/7) * 1.0
        # = 100 / 84 ≈ 1.190476
        expected_rate = 100.0 / 84.0
        @test all(r -> isapprox(r, expected_rate, rtol = 1.0e-6), result.emission_rate)
    end

    @testset "temporal_allocate preserves mass" begin
        # Over a full year with uniform profiles, total hourly rates
        # summed should approximate the annual value when properly scaled
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1.0]
        )
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
            profile_id = [1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), fill(1.0 / 24.0, 24)]
        )
        xref = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        # One full week (168 hours)
        ep_start = DateTime(2019, 7, 1, 0)
        ep_end = DateTime(2019, 7, 8, 0)

        result = temporal_allocate(emissions, profiles, xref, ep_start, ep_end)
        @test nrow(result) == 168  # 7 * 24 hours

        # With uniform profiles, each hour gets ann/84, and a week has 168 hours
        # Total = 168 * (1/84) = 2.0
        total = sum(result.emission_rate)
        @test total ≈ 168.0 / 84.0 rtol = 1.0e-6
    end

    @testset "temporal_allocate empty emissions" begin
        emissions = DataFrame(
            FIPS = String[],
            SCC = String[],
            POLID = String[],
            ANN_VALUE = Float64[]
        )
        profiles = DataFrame(
            profile_type = String[],
            profile_id = Int[],
            factors = Vector{Float64}[]
        )
        xref = DataFrame(
            FIPS = String[],
            SCC = String[],
            monthly_id = Int[],
            weekly_id = Int[],
            diurnal_id = Int[]
        )

        result = temporal_allocate(emissions, profiles, xref, DateTime(2019, 1, 1), DateTime(2019, 1, 2))
        @test nrow(result) == 0
        @test hasproperty(result, :emission_rate)
    end

    @testset "temporal_allocate with Unitful ANN_VALUE" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1.0e-3u"kg/s"]
        )
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
            profile_id = [1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), fill(1.0 / 24.0, 24)]
        )
        xref = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        result = temporal_allocate(
            emissions, profiles, xref,
            DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 1)
        )
        @test nrow(result) == 1
        @test result[1, :emission_rate] ≈ 1.0e-3 / 84.0 rtol = 1.0e-6
    end

    @testset "temporal_allocate multiple sources" begin
        emissions = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0]
        )
        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
            profile_id = [1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), fill(1.0 / 24.0, 24)]
        )
        xref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        result = temporal_allocate(
            emissions, profiles, xref,
            DateTime(2019, 7, 1, 0), DateTime(2019, 7, 1, 2)
        )
        @test nrow(result) == 4  # 2 sources * 2 hours
        nox_rows = filter(r -> r.POLID == "NOX", result)
        voc_rows = filter(r -> r.POLID == "VOC", result)
        @test nrow(nox_rows) == 2
        @test nrow(voc_rows) == 2
    end
end
