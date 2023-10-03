@testset "composite sp." begin
    gas = Gas()
    gas.X = Xair = Dict(
        "N2" => 0.78084,
        "Ar" => 0.009365,
        "Air" => 0.0,
        "H2O" => 0.0,
        "CO2" => 0.000319,
        "O2" => 0.209476)

    gas.T = T = IdealGases.Tstd
    gas.P = P = IdealGases.Pstd

    Air = species_in_spdict("Air")
    composite_air = IdealGases.generate_composite_species(gas.X, "air")

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
    
end