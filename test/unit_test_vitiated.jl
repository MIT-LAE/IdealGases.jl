@testset "vitiated sp." begin
    FAR = 0.5
    Air = IdealGases.DryAir
    O2 = species_in_spdict("O2")
    CH4 = species_in_spdict("CH4")

    # Oxygen combustion
    FARst = IdealGases.stoich_molar_FOR(CH4, O2)
    @test FARst == 0.5

    FAR = 0.2
    massFAR = FAR * CH4.MW / O2.MW
    Xdict = IdealGases.vitiated_mixture(CH4, O2, massFAR)
    @test Xdict["O2"] == (1 - FAR / FARst) / (1 + FAR)

    # Dry Air combustion 
    FARst = IdealGases.stoich_molar_FOR(CH4, Air)
    @test FARst == 0.104738

    FAR = 0.01
    massFAR = FAR * CH4.MW / Air.MW
    Xdict = IdealGases.vitiated_mixture(CH4, Air, massFAR)
    @test Xdict["O2"] == Air.composition["O2"] * (1 - FAR / FARst) / (1 + FAR)

    gas = Gas()
    gas.X = IdealGases.Xair

    ΔX = FAR .* IdealGases.reaction_change_molar_fraction("CH4")
    ΔXdict = Dict(zip(["CO2", "N2", "H2O", "O2"], ΔX))

    Xdict = mergewith(+, gas.Xdict, ΔXdict)
    gas.X = Xdict

    T, P = IdealGases.Tstd, IdealGases.Pstd
    gas.T = T
    gas.P = P
    burntair = IdealGases.vitiated_species(CH4, Air, massFAR)

    #Gas mixture properties
    @test burntair.MW ≈ gas.MW

    #Thermo 
    @test gas.cp ≈ IdealGases.Cp(T, burntair)
    @test gas.h ≈ IdealGases.h(T, burntair)
    @test gas.s ≈ IdealGases.s(T, P, burntair)

    gas.T = 2 * T
    gas.P = 2 * P

    @test gas.cp ≈ IdealGases.Cp(2 * T, burntair)
    @test gas.h ≈ IdealGases.h(2 * T, burntair)
    @test gas.s ≈ IdealGases.s(2 * T, 2 * P, burntair)

    #Test a wide range of T,P
    for Ti in range(100.0, 3000.0, step = 200.0)
        gas.T = Ti
        @test gas.cp ≈ IdealGases.Cp(gas.T, burntair)
        @test gas.h ≈ IdealGases.h(gas.T, burntair)
        for Pi in range(P, 30 * P, step = 5 * P)
            gas.P = Pi
            @test gas.s ≈ IdealGases.s(gas.TP..., burntair)
        end
    end


end

@testset "Gas1D" begin

    #Test with Dry Air first
    gas = Gas()
    gas1 = Gas1D()

    gas.X = DryAir.composition
    @test gas.MW ≈ gas1.comp_sp.MW
    gas.T = gas1.T = 2000.0
    @test gas.T == gas1.T
    @test gas.cp ≈ gas1.cp
    @test gas.h ≈ gas1.h
    @test gas.s ≈ gas1.s

    gas.h = gas1.h = -2000.0
    @test gas.T ≈ gas1.T
    @test gas.cp ≈ gas1.cp
    @test gas.h ≈ gas1.h
    @test gas.s ≈ gas1.s

    Xburnt = IdealGases.vitiated_mixture("CH4", "Air", 0.05)
    burntgas = IdealGases.vitiated_species("CH4", "Air", 0.05)

    gas.X = Xburnt
    gas1.comp_sp = burntgas
    #Test composition first
    @test gas.MW ≈ gas1.comp_sp.MW
    for (key, val) in gas.Xdict
        @test gas.Xdict[key] ≈ gas1.comp_sp.composition[key]
    end

    Tstd, Pstd = IdealGases.Tstd, IdealGases.Pstd
    #Test a wide range of T,P
    for Ti in range(100.0, 3000.0, length = 5)
        gas.T = gas1.T = Ti
        @test gas.T == gas1.T
        @test gas.cp ≈ gas1.cp
        @test gas.h ≈ gas1.h
        for Pi in range(Pstd, 30 * Pstd, length = 5)
            gas.P = gas1.P = Pi
            @test gas.s ≈ gas1.s
        end
    end

    gas.P = gas1.P = Pstd
    for hi in range(-5e3, 5e6, length = 10)
        gas.h = gas1.h = hi
        @test gas.h ≈ gas1.h
        @test gas.T ≈ gas1.T
        @test gas.cp ≈ gas1.cp
        @test gas.s ≈ gas1.s
    end

    #Test setting T,P together
    set_TP!(gas, Tstd, Pstd)
    set_TP!(gas1, Tstd, Pstd)
    @test gas.h ≈ gas1.h
    @test gas.T ≈ gas1.T
    @test gas.cp ≈ gas1.cp
    @test gas.s ≈ gas1.s

    #Test finite Δh increase at a given poly eff
    ηp = 0.95
    set_Δh!(gas, 200.0, ηp)
    set_Δh!(gas1, 200.0, ηp)
    @test gas.h ≈ gas1.h
    @test gas.T ≈ gas1.T
    @test gas.cp ≈ gas1.cp
    @test gas.s ≈ gas1.s

    #Test compression at given poly eff
    set_TP!(gas, Tstd, Pstd)
    set_TP!(gas1, Tstd, Pstd)
    PR = 10.0
    IdealGases.compress(gas, PR, ηp)
    IdealGases.compress(gas1, PR, ηp)
    @test gas.P ≈ Pstd * PR
    @test gas1.P ≈ Pstd * PR
    @test gas.h ≈ gas1.h
    @test gas.T ≈ gas1.T
    @test gas.cp ≈ gas1.cp
    @test gas.s ≈ gas1.s

    # Test expansion
    set_TP!(gas, Tstd, Pstd)
    set_TP!(gas1, Tstd, Pstd)
    PR = 0.5
    IdealGases.expand(gas, PR, ηp)
    IdealGases.expand(gas1, PR, ηp)
    @test gas.P ≈ Pstd * PR
    @test gas1.P ≈ Pstd * PR
    @test gas.T < Tstd
    @test gas1.T < Tstd
    @test gas.h ≈ gas1.h
    @test gas.T ≈ gas1.T
    @test gas.cp ≈ gas1.cp
    @test gas.s ≈ gas1.s
end
