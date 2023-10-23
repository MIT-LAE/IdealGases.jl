@testset "combustion" begin
    CH4 = species_in_spdict("CH4")
    @test IdealGases.LHV(CH4) == 5.002736488044851e7
    @test IdealGases.LHV("CH4") == 5.002736488044851e7

    @test IdealGases.AFT(CH4) == 2200.638532182601#2376.6102617357988
    O2 = species_in_spdict("O2")
    @test IdealGases.AFT(CH4, O2) == 5280.53933225877
end