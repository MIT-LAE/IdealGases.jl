@testset "vitiated sp." begin
    FAR = 0.5
    Air = IdealGases.DryAir
    O2 = species_in_spdict("O2")
    CH4 = species_in_spdict("CH4")

    # Oxygen combustion
    FARst = IdealGases.stoich_molar_FOR(CH4, O2)
    @test FARst == 0.5

    FAR = 0.2
    massFAR = FAR*CH4.MW/O2.MW
    Xdict = IdealGases.vitiated_mixture(CH4, O2, massFAR)
    @test Xdict["O2"] == 1 - FAR/FARst

    # Dry Air combustion 
    FARst = IdealGases.stoich_molar_FOR(CH4, Air)
    @test FARst == 0.104738

    FAR = 0.01
    massFAR = FAR*CH4.MW/Air.MW
    Xdict = IdealGases.vitiated_mixture(CH4, Air, massFAR)
    @test Xdict["O2"] == Air.composition["O2"]*(1 - FAR/FARst)

    gas = Gas()
    gas.X = IdealGases.Xair

    ΔX = FAR.*IdealGases.reaction_change_molar_fraction("CH4")
    ΔXdict = Dict(zip(["CO2", "N2", "H2O", "O2"], ΔX))

    Xdict = mergewith(+, gas.Xdict, ΔXdict)
    gas.X = Xdict

    T,P = IdealGases.Tstd, IdealGases.Pstd
    gas.T = T; gas.P = P
    burntair = IdealGases.vitiated_species(CH4, Air, massFAR)

    #Gas mixture properties
    @test burntair.MW ≈ gas.MW

    #Thermo 
    @test gas.cp ≈ IdealGases.Cp(T, burntair)
    @test gas.h ≈ IdealGases.h(T, burntair)
    @test gas.s ≈ IdealGases.s(T, P, burntair)

    gas.T = 2*T; gas.P = 2*P

    @test gas.cp ≈ IdealGases.Cp(2*T, burntair)
    @test gas.h ≈ IdealGases.h(2*T, burntair)
    @test gas.s ≈ IdealGases.s(2*T, 2*P, burntair)

    #Test a wide range of T,P
    for Ti in range(100.0, 3000.0, step = 200.0)
        gas.T = Ti
        @test gas.cp ≈ IdealGases.Cp(gas.T, burntair)
        @test gas.h ≈ IdealGases.h(gas.T, burntair)
        for Pi in range(P, 30*P, step = 5*P)
            gas.P = Pi
            @test gas.s ≈ IdealGases.s(gas.TP..., burntair)
        end
    end


end