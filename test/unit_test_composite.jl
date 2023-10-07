@testset "composite sp." begin
    gas = Gas()
    gas.X = Xair = IdealGases.Xair

    gas.T = T = IdealGases.Tstd
    gas.P = P = IdealGases.Pstd

    Air = species_in_spdict("Air")
    composite_air = IdealGases.generate_composite_species(gas.X, "air")

    @test IdealGases.Cp(T, composite_air) ≈ gas.cp
    @test IdealGases.h(T, composite_air) ≈ gas.h
    @test IdealGases.s(T, P, composite_air) ≈ gas.s

    @test composite_air.MW ≈ Air.MW atol=1e-7 #that's the no. of decimals in thermo.inp
    @test composite_air.Hf ≈ Air.Hf atol=1e-2 #that's the no. of decimals in thermo.inp
    for i in eachindex(Air.alow)
        @test composite_air.alow[i] ≈ Air.alow[i] rtol=1e-6
        @test composite_air.ahigh[i] ≈ Air.ahigh[i] rtol=1e-4
    end

    @test IdealGases.Cp(T, Air) ≈ IdealGases.Cp(T, composite_air) rtol=1e-8
    @test IdealGases.h(T, Air) ≈ IdealGases.h(T, composite_air) rtol=1e-6
    @test IdealGases.s(T, P, Air) ≈ IdealGases.s(T, P, composite_air) rtol=1e-8
    P *= 2
    T *= 5
    @test IdealGases.s(T, P, Air) ≈ IdealGases.s(T, P, composite_air) rtol=1e-8
    
    for Ti in rand(200.0:2000.0,5)
        gas.T = Ti
        @test gas.cp ≈ IdealGases.Cp(Ti, composite_air)
        @test gas.h ≈ IdealGases.h(Ti, composite_air)
        for Pi in rand(101325.0:5*101325.0,5)
            gas.P = Pi
            @test gas.s≈ IdealGases.s(Ti, gas.P, composite_air)
            @test IdealGases.s(gas.TP..., Air) ≈ IdealGases.s(gas.TP..., composite_air)
        end
    end
end