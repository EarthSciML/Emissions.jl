using DataFrames, Unitful

@testset "Validation tests" begin
    @testset "check_duplicates finds duplicates" begin
        df = DataFrame(
            FIPS = ["36001", "36001", "36005"],
            SCC = ["2103007000", "2103007000", "2103007000"],
            POLID = ["NOX", "NOX", "NOX"],
            ANN_VALUE = [100.0, 50.0, 200.0]
        )
        dups = check_duplicates(df)
        @test nrow(dups) == 1
        @test dups[1, :dup_count] == 2
        @test dups[1, :FIPS] == "36001"
    end

    @testset "check_duplicates no duplicates" begin
        df = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0]
        )
        dups = check_duplicates(df)
        @test nrow(dups) == 0
    end

    @testset "check_duplicates empty DataFrame" begin
        df = DataFrame(FIPS = String[], SCC = String[], POLID = String[], ANN_VALUE = Float64[])
        dups = check_duplicates(df)
        @test nrow(dups) == 0
    end

    @testset "check_ranges negative values" begin
        df = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            ANN_VALUE = [-100.0, 50.0]
        )
        issues = check_ranges(df)
        @test nrow(issues) == 1
        @test issues[1, :range_issue] == :negative
    end

    @testset "check_ranges zero values" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [0.0]
        )
        issues = check_ranges(df)
        @test nrow(issues) == 1
        @test issues[1, :range_issue] == :zero
    end

    @testset "check_ranges extreme values" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [2e12]
        )
        issues = check_ranges(df)
        @test nrow(issues) == 1
        @test issues[1, :range_issue] == :extreme
    end

    @testset "check_ranges custom max_value" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [500.0]
        )
        issues = check_ranges(df; max_value = 100.0)
        @test nrow(issues) == 1
        @test issues[1, :range_issue] == :extreme
    end

    @testset "check_ranges all valid" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0]
        )
        issues = check_ranges(df)
        @test nrow(issues) == 0
    end

    @testset "check_ranges with Unitful values" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [-1.0e-3u"kg/s"]
        )
        issues = check_ranges(df)
        @test nrow(issues) == 1
        @test issues[1, :range_issue] == :negative
    end

    @testset "check_completeness missing fields" begin
        df = DataFrame(
            FIPS = ["36001", missing],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            ANN_VALUE = [100.0, 50.0]
        )
        issues = check_completeness(df)
        @test nrow(issues) == 1
        @test :FIPS in issues[1, :missing_fields]
    end

    @testset "check_completeness empty strings" begin
        df = DataFrame(
            FIPS = ["36001", ""],
            SCC = ["2103007000", "2103007000"],
            POLID = ["NOX", "NOX"],
            ANN_VALUE = [100.0, 50.0]
        )
        issues = check_completeness(df)
        @test nrow(issues) == 1
        @test :FIPS in issues[1, :missing_fields]
    end

    @testset "check_completeness all complete" begin
        df = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0]
        )
        issues = check_completeness(df)
        @test nrow(issues) == 0
    end

    @testset "validate_inventory clean data" begin
        df = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007001"],
            POLID = ["NOX", "VOC"],
            ANN_VALUE = [100.0, 50.0]
        )
        result = validate_inventory(df)
        @test result.valid == true
        @test isempty(result.errors)
        @test result.n_duplicates == 0
        @test result.n_range_issues == 0
        @test result.n_missing_fields == 0
    end

    @testset "validate_inventory with issues" begin
        df = DataFrame(
            FIPS = ["36001", "36001", missing],
            SCC = ["2103007000", "2103007000", "2103007001"],
            POLID = ["NOX", "NOX", "VOC"],
            ANN_VALUE = [-100.0, 50.0, 200.0]
        )
        result = validate_inventory(df)
        @test result.valid == false
        @test result.n_duplicates == 1
        @test result.n_range_issues == 1
        @test result.n_missing_fields == 1
        @test !isempty(result.errors)
        @test !isempty(result.warnings)
    end

    @testset "validate_inventory empty DataFrame" begin
        df = DataFrame(FIPS = String[], SCC = String[], POLID = String[], ANN_VALUE = Float64[])
        result = validate_inventory(df)
        @test result.valid == true
        @test result.n_duplicates == 0
    end
end
