using SparseArrays

@testset "Surrogates tests" begin
    @testset "generate_countySurrogate" begin
        data = Dict(
            "01001" => sparse([1, 2], [1, 1], [0.5, 0.5], 2, 2),
            "01002" => sparse([1], [2], [1.0], 2, 2),
        )
        weight = sparse([1, 2], [1, 2], [2.0, 3.0], 2, 2)
        grid = sparse([1, 1, 2, 2], [1, 2, 1, 2], [1.0, 1.0, 1.0, 1.0], 2, 2)

        result = generate_countySurrogate(data, weight, grid)

        @test haskey(result, "01001")
        @test haskey(result, "01002")

        # Check that fractions sum to 1 for county 01001
        @test sum(result["01001"]) ≈ 1.0
    end

    @testset "generate_countySurrogate zero weight fallback" begin
        data = Dict(
            "01001" => sparse([1], [1], [1.0], 2, 2),
        )
        weight = spzeros(2, 2)
        grid = sparse([1, 1, 2, 2], [1, 2, 1, 2], [1.0, 1.0, 1.0, 1.0], 2, 2)

        result = generate_countySurrogate(data, weight, grid)
        @test sum(result["01001"]) ≈ 1.0
    end

    @testset "update_locIndex with matching FIPS" begin
        locIndex = Dict{String,IndexInfo}()
        surrogates = Dict(
            "01001" => sparse([1, 2], [1, 2], [0.6, 0.4], 3, 3),
        )

        idx = update_locIndex(locIndex, "01001", surrogates)
        @test idx.inGrid == true
        @test length(idx.rows) == 2
        @test haskey(locIndex, "01001")
    end

    @testset "update_locIndex without matching FIPS" begin
        locIndex = Dict{String,IndexInfo}()
        surrogates = Dict{String,SparseMatrixCSC{Float64,Int}}()

        idx = update_locIndex(locIndex, "99999", surrogates)
        @test idx.inGrid == false
        @test isempty(idx.rows)
    end

    @testset "find_column_name" begin
        using DataFrames
        df = DataFrame(Population=[1, 2], Area=[3.0, 4.0])
        @test Emissions.find_column_name(df, "population") == "Population"
        @test Emissions.find_column_name(df, "POPULATION") == "Population"
        @test_throws ErrorException Emissions.find_column_name(df, "nonexistent")
    end
end
