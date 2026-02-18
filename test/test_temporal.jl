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
        # With uniform profiles: rate = ann * (mf*12) * wf * (df*24)
        # With uniform: mf=1/12, wf=1.0, df=1/24
        # rate = 100 * (1/12 * 12) * 1.0 * (1/24 * 24) = 100 * 1 * 1 * 1 = 100.0
        expected_rate = 100.0
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

        # With uniform profiles, each hour gets ann * 1.0, and a week has 168 hours
        # Total = 168 * 1.0 = 168.0
        total = sum(result.emission_rate)
        @test total ≈ 168.0 rtol = 1.0e-6
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
        # With uniform profiles, hourly rate = annual rate
        @test result[1, :emission_rate] ≈ 1.0e-3 rtol = 1.0e-6
    end

    @testset "temporal_allocate per-FIPS timezone" begin
        # Two sources in different timezones: FIPS 36001 (Eastern, UTC-5)
        # and FIPS 06001 (Pacific, UTC-8)
        emissions = DataFrame(
            FIPS = ["36001", "06001"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            ANN_VALUE = [100.0, 100.0]
        )
        # Non-uniform diurnal profile: daytime hours (8-17) get 2x weight
        diurnal = zeros(24)
        for h in 1:24
            diurnal[h] = (h >= 9 && h <= 18) ? 2.0 / 24.0 : 0.4 / 24.0
        end
        # Normalize so sum = 1.0
        diurnal ./= sum(diurnal)

        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "DIURNAL"],
            profile_id = [1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), diurnal]
        )
        xref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        # Test at UTC hour 14 (9am Eastern, 6am Pacific)
        ep_start = DateTime(2019, 7, 1, 14)
        ep_end = DateTime(2019, 7, 1, 15)

        timezone_map = Dict("36001" => -5, "06001" => -8)
        result = temporal_allocate(emissions, profiles, xref, ep_start, ep_end;
                                   timezone_map = timezone_map)
        @test nrow(result) == 2

        # Eastern source: local hour 9 (daytime), Pacific source: local hour 6 (night)
        eastern = filter(r -> r.FIPS == "36001", result)
        pacific = filter(r -> r.FIPS == "06001", result)
        @test nrow(eastern) == 1
        @test nrow(pacific) == 1
        # The eastern source should have higher emissions (daytime profile)
        @test eastern[1, :emission_rate] > pacific[1, :emission_rate]
    end

    @testset "temporal_allocate day-specific diurnal profiles" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0]
        )
        # Weekday profile: concentrated in work hours
        weekday_diurnal = fill(0.01, 24)
        weekday_diurnal[9:17] .= 0.09  # 9 hours * 0.09 + 15 hours * 0.01 = 0.81 + 0.15 = 0.96
        weekday_diurnal ./= sum(weekday_diurnal)

        # Weekend profile: uniform
        weekend_diurnal = fill(1.0 / 24.0, 24)

        profiles = DataFrame(
            profile_type = ["MONTHLY", "WEEKLY", "WEEKDAY", "WEEKEND"],
            profile_id = [1, 1, 1, 1],
            factors = [fill(1.0 / 12.0, 12), fill(1.0, 7), weekday_diurnal, weekend_diurnal]
        )
        xref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        # Monday July 1, 2019 at hour 12 (work hours) - should use WEEKDAY profile
        result_weekday = temporal_allocate(emissions, profiles, xref,
            DateTime(2019, 7, 1, 12), DateTime(2019, 7, 1, 13))
        @test nrow(result_weekday) == 1

        # Saturday July 6, 2019 at hour 12 - should use WEEKEND profile
        result_weekend = temporal_allocate(emissions, profiles, xref,
            DateTime(2019, 7, 6, 12), DateTime(2019, 7, 6, 13))
        @test nrow(result_weekend) == 1

        # Weekday hour 12 should differ from weekend hour 12
        @test result_weekday[1, :emission_rate] != result_weekend[1, :emission_rate]
    end

    @testset "temporal_allocate full year mass conservation" begin
        # With uniform profiles, summing hourly rates over a full year
        # and dividing by 8760 should recover the annual value
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [1000.0]
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

        # Full year
        ep_start = DateTime(2019, 1, 1, 0)
        ep_end = DateTime(2020, 1, 1, 0)
        n_hours = 8760

        result = temporal_allocate(emissions, profiles, xref, ep_start, ep_end)
        @test nrow(result) == n_hours

        # With uniform profiles: each hour gets ann * (1/12*12) * 1.0 * (1/24*24) = ann
        # The per-hour rate should be constant and equal to the annual rate.
        expected_rate = 1000.0
        @test all(r -> isapprox(r, expected_rate, rtol = 1e-6), result.emission_rate)
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

    @testset "read_day_specific" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Day-specific emissions")
            println(io, "036001;2103007000;NOX;07/01/2019;500.0")
            println(io, "036001;2103007000;NOX;07/02/2019;600.0!daily override")
        end
        ds = read_day_specific(tmpfile)
        @test nrow(ds) == 2
        @test ds[1, :FIPS] == "36001"
        @test ds[1, :SCC] == "2103007000"
        @test ds[1, :POLID] == "NOX"
        @test ds[1, :date] == Date(2019, 7, 1)
        @test ds[1, :day_value] ≈ 500.0
        @test ds[2, :day_value] ≈ 600.0
        rm(tmpfile)
    end

    @testset "read_hour_specific" begin
        tmpfile = tempname()
        hourly_vals = join(fill("10.0", 24), ";")
        open(tmpfile, "w") do io
            println(io, "# Hour-specific emissions")
            println(io, "036001;2103007000;NOX;07/01/2019;$hourly_vals")
        end
        hs = read_hour_specific(tmpfile)
        @test nrow(hs) == 1
        @test hs[1, :FIPS] == "36001"
        @test hs[1, :date] == Date(2019, 7, 1)
        @test length(hs[1, :hourly_values]) == 24
        @test all(v -> v ≈ 10.0, hs[1, :hourly_values])
        rm(tmpfile)
    end

    @testset "temporal_allocate with hour_specific override" begin
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
            FIPS = ["00000"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        # Hour-specific: override with known values
        hour_specific = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            date = [Date(2019, 7, 1)],
            hourly_values = [collect(1.0:24.0)]
        )

        ep_start = DateTime(2019, 7, 1, 0)
        ep_end = DateTime(2019, 7, 1, 3)  # 3 hours

        result = temporal_allocate(emissions, profiles, xref, ep_start, ep_end;
            hour_specific = hour_specific)
        @test nrow(result) == 3
        # Hour 0 -> hour_idx 1 -> value 1.0
        # Hour 1 -> hour_idx 2 -> value 2.0
        # Hour 2 -> hour_idx 3 -> value 3.0
        sort!(result, :hour)
        @test result[1, :emission_rate] ≈ 1.0
        @test result[2, :emission_rate] ≈ 2.0
        @test result[3, :emission_rate] ≈ 3.0
    end

    @testset "temporal_allocate with day_specific override" begin
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
            FIPS = ["00000"],
            SCC = ["2103007000"],
            monthly_id = [1],
            weekly_id = [1],
            diurnal_id = [1]
        )

        # Day-specific: override daily total
        day_specific = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            date = [Date(2019, 7, 1)],
            day_value = [240.0]
        )

        ep_start = DateTime(2019, 7, 1, 0)
        ep_end = DateTime(2019, 7, 1, 3)

        result = temporal_allocate(emissions, profiles, xref, ep_start, ep_end;
            day_specific = day_specific)
        @test nrow(result) == 3
        # With uniform diurnal (1/24), day_value * (1/24) * 24 = 240.0 * 1.0 = 240.0
        # Wait: hourly_rate = day_val * df * 24 = 240 * (1/24) * 24 = 240.0
        @test all(r -> r ≈ 240.0, result.emission_rate)
    end

    @testset "temporal_allocate backward compatibility" begin
        # Verify that calling without day/hour specific produces same results
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
        ep_end = DateTime(2019, 7, 1, 3)

        result_old = temporal_allocate(emissions, profiles, xref, ep_start, ep_end)
        result_new = temporal_allocate(emissions, profiles, xref, ep_start, ep_end;
            day_specific = DataFrame(), hour_specific = DataFrame())

        @test nrow(result_old) == nrow(result_new)
        @test all(result_old.emission_rate .≈ result_new.emission_rate)
    end
end
