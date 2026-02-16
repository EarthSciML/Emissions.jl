using DataFrames

@testset "Biogenic tests" begin
    @testset "read_beld" begin
        grid = NewGridIrregular("test", 3, 3, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# BELD land use data")
            println(io, "1,Deciduous Forest,0.5")
            println(io, "1,Cropland,0.3")
            println(io, "5,Coniferous Forest,0.8")
        end
        beld = read_beld(tmpfile, grid)
        @test nrow(beld) == 3
        @test beld[1, :cell_index] == 1
        @test beld[1, :land_use_type] == "Deciduous Forest"
        @test beld[1, :fraction] ≈ 0.5
        @test beld[3, :cell_index] == 5
        rm(tmpfile)
    end

    @testset "read_beld skips out of range" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "1,Forest,0.5")
            println(io, "10,Forest,0.5")  # Out of range for 2x2 grid
            println(io, "0,Forest,0.5")   # Out of range
        end
        beld = read_beld(tmpfile, grid)
        @test nrow(beld) == 1
        rm(tmpfile)
    end

    @testset "read_emission_factors" begin
        tmpfile = tempname()
        open(tmpfile, "w") do io
            println(io, "# Emission factors")
            println(io, "Deciduous Forest,ISOP,10000.0,1000.0")
            println(io, "Deciduous Forest,TERP,500.0,100.0")
            println(io, "Coniferous Forest,ISOP,5000.0,500.0")
            println(io, "Coniferous Forest,TERP,2000.0,400.0")
        end
        ef = read_emission_factors(tmpfile)
        @test nrow(ef) == 4
        @test ef[1, :land_use_type] == "Deciduous Forest"
        @test ef[1, :species] == "ISOP"
        @test ef[1, :summer_factor] ≈ 10000.0
        @test ef[1, :winter_factor] ≈ 1000.0
        rm(tmpfile)
    end

    @testset "temperature_adjustment isoprene at standard conditions" begin
        # At 30°C (303.15K), the adjustment should be close to reference
        γ = temperature_adjustment(303.15, "ISOP")
        @test γ > 0.0
        @test isfinite(γ)
    end

    @testset "temperature_adjustment isoprene increases with temperature" begin
        γ_low = temperature_adjustment(290.0, "ISOP")
        γ_mid = temperature_adjustment(303.15, "ISOP")
        γ_high = temperature_adjustment(310.0, "ISOP")
        # Should increase from low to mid temperature
        @test γ_mid > γ_low
    end

    @testset "temperature_adjustment isoprene decreases at extreme heat" begin
        # Isoprene has a bell curve peaking around 314K
        γ_opt = temperature_adjustment(314.0, "ISOP")
        γ_hot = temperature_adjustment(340.0, "ISOP")
        @test γ_opt > γ_hot
    end

    @testset "temperature_adjustment terpene at standard conditions" begin
        γ = temperature_adjustment(303.15, "TERP")
        @test γ ≈ 1.0
    end

    @testset "temperature_adjustment terpene increases with temperature" begin
        γ_low = temperature_adjustment(280.0, "TERP")
        γ_high = temperature_adjustment(310.0, "TERP")
        @test γ_high > γ_low
    end

    @testset "temperature_adjustment cold temperature" begin
        # Very cold should give very low emission
        γ_isop = temperature_adjustment(250.0, "ISOP")
        γ_terp = temperature_adjustment(250.0, "TERP")
        @test γ_isop >= 0.0
        @test γ_terp > 0.0
        @test γ_terp < 1.0
    end

    @testset "light_adjustment reference conditions" begin
        # At PAR=1000 (standard), should be close to reference
        γ = light_adjustment(1000.0)
        @test γ > 0.0
        @test isfinite(γ)
    end

    @testset "light_adjustment zero PAR" begin
        @test light_adjustment(0.0) == 0.0
        @test light_adjustment(-10.0) == 0.0
    end

    @testset "light_adjustment increases with PAR" begin
        γ_low = light_adjustment(200.0)
        γ_high = light_adjustment(1000.0)
        @test γ_high > γ_low
    end

    @testset "light_adjustment saturates" begin
        γ_high = light_adjustment(1000.0)
        γ_very_high = light_adjustment(2000.0)
        # Increase should be diminishing
        @test γ_very_high > γ_high
        @test (γ_very_high - γ_high) < (γ_high - light_adjustment(0.0))
    end

    @testset "compute_biogenic_emissions basic" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        beld_file = tempname()
        open(beld_file, "w") do io
            println(io, "1,Deciduous Forest,0.5")
            println(io, "2,Coniferous Forest,0.3")
        end

        ef_file = tempname()
        open(ef_file, "w") do io
            println(io, "Deciduous Forest,ISOP,10000.0,1000.0")
            println(io, "Deciduous Forest,TERP,500.0,100.0")
            println(io, "Coniferous Forest,TERP,2000.0,400.0")
        end

        config = BiogenicConfig(beld_file, ef_file, :summer)
        temperature = fill(303.15, 4)
        par = fill(1000.0, 4)

        result = compute_biogenic_emissions(config, grid, temperature, par)
        @test nrow(result) > 0
        @test hasproperty(result, :grid_row)
        @test hasproperty(result, :grid_col)
        @test hasproperty(result, :species)
        @test hasproperty(result, :emission_rate)
        @test all(result.emission_rate .> 0.0)

        rm(beld_file)
        rm(ef_file)
    end

    @testset "compute_biogenic_emissions correct species" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        beld_file = tempname()
        open(beld_file, "w") do io
            println(io, "1,Forest,1.0")
        end

        ef_file = tempname()
        open(ef_file, "w") do io
            println(io, "Forest,ISOP,10000.0,1000.0")
            println(io, "Forest,TERP,500.0,100.0")
        end

        config = BiogenicConfig(beld_file, ef_file, :summer)
        temperature = fill(303.15, 4)
        par = fill(1000.0, 4)

        result = compute_biogenic_emissions(config, grid, temperature, par)
        species = unique(result.species)
        @test "ISOP" in species
        @test "TERP" in species

        rm(beld_file)
        rm(ef_file)
    end

    @testset "compute_biogenic_emissions winter season" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        beld_file = tempname()
        open(beld_file, "w") do io
            println(io, "1,Forest,1.0")
        end

        ef_file = tempname()
        open(ef_file, "w") do io
            println(io, "Forest,TERP,500.0,100.0")
        end

        config_summer = BiogenicConfig(beld_file, ef_file, :summer)
        config_winter = BiogenicConfig(beld_file, ef_file, :winter)
        temperature = fill(303.15, 4)
        par = fill(1000.0, 4)

        result_summer = compute_biogenic_emissions(config_summer, grid, temperature, par)
        result_winter = compute_biogenic_emissions(config_winter, grid, temperature, par)

        # Winter should have lower emissions (factor 100 vs 500)
        @test sum(result_winter.emission_rate) < sum(result_summer.emission_rate)

        rm(beld_file)
        rm(ef_file)
    end

    @testset "compute_biogenic_emissions mismatched grid size errors" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        beld_file = tempname()
        open(beld_file, "w") do io
            println(io, "1,Forest,1.0")
        end
        ef_file = tempname()
        open(ef_file, "w") do io
            println(io, "Forest,ISOP,10000.0,1000.0")
        end

        config = BiogenicConfig(beld_file, ef_file, :summer)

        # Wrong length temperature vector
        @test_throws ArgumentError compute_biogenic_emissions(
            config, grid, fill(300.0, 3), fill(1000.0, 4)
        )

        rm(beld_file)
        rm(ef_file)
    end

    @testset "compute_biogenic_emissions empty land use" begin
        grid = NewGridIrregular("test", 2, 2, "EPSG:4326", 1.0, 1.0, 0.0, 0.0)

        beld_file = tempname()
        open(beld_file, "w") do io
            println(io, "# Empty")
        end
        ef_file = tempname()
        open(ef_file, "w") do io
            println(io, "Forest,ISOP,10000.0,1000.0")
        end

        config = BiogenicConfig(beld_file, ef_file, :summer)
        result = compute_biogenic_emissions(config, grid, fill(300.0, 4), fill(1000.0, 4))
        @test nrow(result) == 0

        rm(beld_file)
        rm(ef_file)
    end
end
