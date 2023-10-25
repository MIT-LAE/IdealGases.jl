@testset "combustion" begin
    CH4 = species_in_spdict("CH4")
    @test IdealGases.LHV(CH4) == 5.002736488044851e7
    @test IdealGases.LHV("CH4") == 5.002736488044851e7

    @test IdealGases.AFT(CH4) == 2376.6102617357988
    O2 = species_in_spdict("O2")
    @test IdealGases.AFT(CH4, O2) == 5280.53933225877
end
@testset "burnt gas funcs." begin
    CH4 = species_in_spdict("CH4")

    FAR = 0.05
    gas1 = IdealGases.vitiated_species(CH4, DryAir, FAR)

    burntgas = IdealGases.fixed_fuel_vitiated_species(CH4, DryAir)
    gas2 = burntgas(FAR)

    @test gas1.MW ≈ gas2.MW
    for (key,val) in gas1.composition
        @test gas1.composition[key] ≈ gas2.composition[key]
    end

end