@testset "composite sp." begin
    gas = Gas()
    gas.X = Xair = IdealGasThermo.Xair

    gas.T = T = IdealGasThermo.Tstd
    gas.P = P = IdealGasThermo.Pstd

    Air = species_in_spdict("Air")
    composite_air = IdealGasThermo.generate_composite_species(gas.X, "air")

    @test IdealGasThermo.Cp(T, composite_air) ≈ gas.cp
    @test IdealGasThermo.h(T, composite_air) ≈ gas.h
    @test IdealGasThermo.s(T, P, composite_air) ≈ gas.s

    @test composite_air.MW ≈ Air.MW atol = 1e-7 #that's the no. of decimals in thermo.inp
    @test composite_air.Hf ≈ Air.Hf atol = 1e-2 #that's the no. of decimals in thermo.inp
    for i in eachindex(Air.alow)
        @test composite_air.alow[i] ≈ Air.alow[i] rtol = 1e-6
        @test composite_air.ahigh[i] ≈ Air.ahigh[i] rtol = 1e-4
    end

    @test IdealGasThermo.Cp(T, Air) ≈ IdealGasThermo.Cp(T, composite_air) rtol = 1e-8
    @test IdealGasThermo.h(T, Air) ≈ IdealGasThermo.h(T, composite_air) rtol = 1e-6
    @test IdealGasThermo.s(T, P, Air) ≈ IdealGasThermo.s(T, P, composite_air) rtol = 1e-8
    P *= 2
    T *= 5
    @test IdealGasThermo.s(T, P, Air) ≈ IdealGasThermo.s(T, P, composite_air) rtol = 1e-8

    for Ti in rand(200.0:2000.0, 5)
        gas.T = Ti
        @test gas.cp ≈ IdealGasThermo.Cp(Ti, composite_air)
        @test gas.h ≈ IdealGasThermo.h(Ti, composite_air)
        for Pi in rand(101325.0:5*101325.0, 5)
            gas.P = Pi
            @test gas.s ≈ IdealGasThermo.s(Ti, gas.P, composite_air)
            @test IdealGasThermo.s(gas.TP..., Air) ≈ IdealGasThermo.s(gas.TP..., composite_air)
        end
    end
end
