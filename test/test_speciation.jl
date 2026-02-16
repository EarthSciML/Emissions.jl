using DataFrames

@testset "Speciation tests" begin
    @testset "read_gspro" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# GSPRO speciation profile file")
            println(io, "P001;VOC;FORM;0.5;1.0;0.03")
            println(io, "P001;VOC;ALD2;0.3;1.0;0.02")
            println(io, "P001;VOC;PAR;0.2;1.0;0.95")
            println(io, "P002;NOX;NO;0.9;1.0;0.90 !inline comment")
            println(io, "P002;NOX;NO2;0.1;1.0;0.10")
        end
        gspro = read_gspro(tmpfile)
        @test nrow(gspro) == 5
        @test gspro[1, :profile_code] == "P001"
        @test gspro[1, :pollutant_id] == "VOC"
        @test gspro[1, :species_id] == "FORM"
        @test gspro[1, :split_factor] ≈ 0.5
        @test gspro[1, :divisor] ≈ 1.0
        @test gspro[1, :mass_fraction] ≈ 0.03
        @test gspro[4, :species_id] == "NO"
        @test gspro[4, :mass_fraction] ≈ 0.9
        rm(tmpfile)
    end

    @testset "read_gspro empty file" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Only comments")
        end
        gspro = read_gspro(tmpfile)
        @test nrow(gspro) == 0
        @test hasproperty(gspro, :profile_code)
        rm(tmpfile)
    end

    @testset "read_gsref" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# GSREF cross-reference file")
            println(io, "036001;2103007000;VOC;P001")
            println(io, "000000;2103007000;NOX;P002 !national default")
        end
        gsref = read_gsref(tmpfile)
        @test nrow(gsref) == 2
        @test gsref[1, :FIPS] == "36001"
        @test gsref[1, :SCC] == "2103007000"
        @test gsref[1, :pollutant_id] == "VOC"
        @test gsref[1, :profile_code] == "P001"
        @test gsref[2, :FIPS] == "00000"
        @test gsref[2, :pollutant_id] == "NOX"
        rm(tmpfile)
    end

    @testset "read_gsref with 5-digit FIPS and short SCC" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "36001;21030;VOC;P001")
        end
        gsref = read_gsref(tmpfile)
        @test nrow(gsref) == 1
        @test gsref[1, :FIPS] == "36001"
        @test gsref[1, :SCC] == "0000021030"
        rm(tmpfile)
    end

    @testset "_match_speciation_profile hierarchical matching" begin
        gsref = DataFrame(
            FIPS = ["36001", "00000", "00000"],
            SCC = ["2103007000", "2103007000", "0000000000"],
            pollutant_id = ["VOC", "VOC", "NOX"],
            profile_code = ["P001", "P002", "P003"],
        )

        # Level 1: Exact FIPS + SCC + pollutant
        result = Emissions._match_speciation_profile(gsref, "36001", "2103007000", "VOC")
        @test result == "P001"

        # Level 2: National default (FIPS=00000) + SCC + pollutant
        result = Emissions._match_speciation_profile(gsref, "99999", "2103007000", "VOC")
        @test result == "P002"

        # Level 3: Pollutant-only (SCC=0000000000)
        result = Emissions._match_speciation_profile(gsref, "99999", "9999999999", "NOX")
        @test result == "P003"

        # No match
        result = Emissions._match_speciation_profile(gsref, "99999", "9999999999", "SO2")
        @test result === nothing
    end

    @testset "build_speciation_matrix mass basis" begin
        emissions = DataFrame(
            FIPS = ["36001", "36001"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["VOC", "NOX"],
            ANN_VALUE = [100.0, 50.0],
        )
        gspro = DataFrame(
            profile_code = ["P001", "P001", "P002", "P002"],
            pollutant_id = ["VOC", "VOC", "NOX", "NOX"],
            species_id = ["FORM", "ALD2", "NO", "NO2"],
            split_factor = [0.5, 0.3, 0.9, 0.1],
            divisor = [1.0, 1.0, 1.0, 1.0],
            mass_fraction = [0.03, 0.02, 0.9, 0.1],
        )
        gsref = DataFrame(
            FIPS = ["00000", "00000"],
            SCC = ["2103007000", "2103007000"],
            pollutant_id = ["VOC", "NOX"],
            profile_code = ["P001", "P002"],
        )

        matrix, species, indices = build_speciation_matrix(emissions, gspro, gsref; basis = :mass)
        @test size(matrix, 1) == 2  # 2 sources
        @test length(species) == 4  # ALD2, FORM, NO, NO2
        @test length(indices) == 2

        # Check VOC row has FORM and ALD2 fractions
        voc_row = indices[1] == 1 ? 1 : 2
        form_col = findfirst(==("FORM"), species)
        ald2_col = findfirst(==("ALD2"), species)
        @test form_col !== nothing
        @test ald2_col !== nothing
    end

    @testset "build_speciation_matrix mole basis" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [50.0],
        )
        gspro = DataFrame(
            profile_code = ["P002", "P002"],
            pollutant_id = ["NOX", "NOX"],
            species_id = ["NO", "NO2"],
            split_factor = [0.9, 0.1],
            divisor = [2.0, 2.0],
            mass_fraction = [0.9, 0.1],
        )
        gsref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            pollutant_id = ["NOX"],
            profile_code = ["P002"],
        )

        matrix, species, indices = build_speciation_matrix(emissions, gspro, gsref; basis = :mole)
        @test size(matrix, 1) == 1
        no_col = findfirst(==("NO"), species)
        @test no_col !== nothing
        @test matrix[1, findfirst(==("NO"), species)] ≈ 0.9 / 2.0
        @test matrix[1, findfirst(==("NO2"), species)] ≈ 0.1 / 2.0
    end

    @testset "speciate_emissions mass basis" begin
        emissions = DataFrame(
            FIPS = ["36001", "36001"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["VOC", "NOX"],
            ANN_VALUE = [100.0, 50.0],
        )
        gspro = DataFrame(
            profile_code = ["P001", "P001", "P002", "P002"],
            pollutant_id = ["VOC", "VOC", "NOX", "NOX"],
            species_id = ["FORM", "ALD2", "NO", "NO2"],
            split_factor = [0.5, 0.3, 0.9, 0.1],
            divisor = [1.0, 1.0, 1.0, 1.0],
            mass_fraction = [0.03, 0.02, 0.9, 0.1],
        )
        gsref = DataFrame(
            FIPS = ["00000", "00000"],
            SCC = ["2103007000", "2103007000"],
            pollutant_id = ["VOC", "NOX"],
            profile_code = ["P001", "P002"],
        )

        result = speciate_emissions(emissions, gspro, gsref; basis = :mass)
        @test hasproperty(result, :species)
        @test !hasproperty(result, :POLID)
        @test nrow(result) == 4  # 2 species for VOC + 2 species for NOX

        # Check VOC -> FORM
        form_rows = filter(r -> r.species == "FORM", result)
        @test nrow(form_rows) == 1
        @test form_rows[1, :ANN_VALUE] ≈ 100.0 * 0.03

        # Check VOC -> ALD2
        ald2_rows = filter(r -> r.species == "ALD2", result)
        @test nrow(ald2_rows) == 1
        @test ald2_rows[1, :ANN_VALUE] ≈ 100.0 * 0.02

        # Check NOX -> NO
        no_rows = filter(r -> r.species == "NO", result)
        @test nrow(no_rows) == 1
        @test no_rows[1, :ANN_VALUE] ≈ 50.0 * 0.9

        # Check NOX -> NO2
        no2_rows = filter(r -> r.species == "NO2", result)
        @test nrow(no2_rows) == 1
        @test no2_rows[1, :ANN_VALUE] ≈ 50.0 * 0.1
    end

    @testset "speciate_emissions mole basis" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0],
        )
        gspro = DataFrame(
            profile_code = ["P002", "P002"],
            pollutant_id = ["NOX", "NOX"],
            species_id = ["NO", "NO2"],
            split_factor = [0.9, 0.1],
            divisor = [2.0, 2.0],
            mass_fraction = [0.9, 0.1],
        )
        gsref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            pollutant_id = ["NOX"],
            profile_code = ["P002"],
        )

        result = speciate_emissions(emissions, gspro, gsref; basis = :mole)
        no_rows = filter(r -> r.species == "NO", result)
        @test no_rows[1, :ANN_VALUE] ≈ 100.0 * 0.9 / 2.0
    end

    @testset "speciate_emissions preserves extra columns" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [50.0],
            COUNTRY = ["USA"],
            LONGITUDE = [-73.5],
            LATITUDE = [40.5],
        )
        gspro = DataFrame(
            profile_code = ["P002"],
            pollutant_id = ["NOX"],
            species_id = ["NO"],
            split_factor = [1.0],
            divisor = [1.0],
            mass_fraction = [0.9],
        )
        gsref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            pollutant_id = ["NOX"],
            profile_code = ["P002"],
        )

        result = speciate_emissions(emissions, gspro, gsref)
        @test hasproperty(result, :COUNTRY)
        @test hasproperty(result, :LONGITUDE)
        @test hasproperty(result, :LATITUDE)
        @test result[1, :COUNTRY] == "USA"
        @test result[1, :LONGITUDE] ≈ -73.5
        @test result[1, :LATITUDE] ≈ 40.5
    end

    @testset "speciate_emissions unmatched records pass through" begin
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["9999999999"],
            POLID = ["CO"],
            ANN_VALUE = [25.0],
        )
        gspro = DataFrame(
            profile_code = String[],
            pollutant_id = String[],
            species_id = String[],
            split_factor = Float64[],
            divisor = Float64[],
            mass_fraction = Float64[],
        )
        gsref = DataFrame(
            FIPS = String[],
            SCC = String[],
            pollutant_id = String[],
            profile_code = String[],
        )

        result = speciate_emissions(emissions, gspro, gsref)
        @test nrow(result) == 1
        @test result[1, :species] == "CO"
        @test result[1, :ANN_VALUE] ≈ 25.0
    end

    @testset "speciate_emissions empty input" begin
        emissions = DataFrame(
            FIPS = String[],
            SCC = String[],
            POLID = String[],
            ANN_VALUE = Float64[],
        )
        gspro = DataFrame(
            profile_code = String[],
            pollutant_id = String[],
            species_id = String[],
            split_factor = Float64[],
            divisor = Float64[],
            mass_fraction = Float64[],
        )
        gsref = DataFrame(
            FIPS = String[],
            SCC = String[],
            pollutant_id = String[],
            profile_code = String[],
        )

        result = speciate_emissions(emissions, gspro, gsref)
        @test nrow(result) == 0
        @test hasproperty(result, :species)
    end

    @testset "speciate_emissions mass conservation" begin
        # NOX with fractions summing to 1.0 should conserve mass
        emissions = DataFrame(
            FIPS = ["36001"],
            SCC = ["2103007000"],
            POLID = ["NOX"],
            ANN_VALUE = [100.0],
        )
        gspro = DataFrame(
            profile_code = ["P002", "P002"],
            pollutant_id = ["NOX", "NOX"],
            species_id = ["NO", "NO2"],
            split_factor = [0.9, 0.1],
            divisor = [1.0, 1.0],
            mass_fraction = [0.9, 0.1],
        )
        gsref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            pollutant_id = ["NOX"],
            profile_code = ["P002"],
        )

        result = speciate_emissions(emissions, gspro, gsref; basis = :mass)
        total = sum(result.ANN_VALUE)
        @test total ≈ 100.0 * (0.9 + 0.1)
    end

    @testset "speciate_emissions multiple sources same pollutant" begin
        emissions = DataFrame(
            FIPS = ["36001", "36005"],
            SCC = ["2103007000", "2103007000"],
            POLID = ["VOC", "VOC"],
            ANN_VALUE = [100.0, 200.0],
        )
        gspro = DataFrame(
            profile_code = ["P001", "P001"],
            pollutant_id = ["VOC", "VOC"],
            species_id = ["FORM", "ALD2"],
            split_factor = [0.5, 0.5],
            divisor = [1.0, 1.0],
            mass_fraction = [0.6, 0.4],
        )
        gsref = DataFrame(
            FIPS = ["00000"],
            SCC = ["2103007000"],
            pollutant_id = ["VOC"],
            profile_code = ["P001"],
        )

        result = speciate_emissions(emissions, gspro, gsref; basis = :mass)
        @test nrow(result) == 4  # 2 sources * 2 species

        # Check county 36001 FORM
        form_36001 = filter(r -> r.species == "FORM" && r.FIPS == "36001", result)
        @test nrow(form_36001) == 1
        @test form_36001[1, :ANN_VALUE] ≈ 100.0 * 0.6

        # Check county 36005 ALD2
        ald2_36005 = filter(r -> r.species == "ALD2" && r.FIPS == "36005", result)
        @test nrow(ald2_36005) == 1
        @test ald2_36005[1, :ANN_VALUE] ≈ 200.0 * 0.4
    end
end
