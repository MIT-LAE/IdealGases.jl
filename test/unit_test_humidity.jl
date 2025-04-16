@testset "humidity" begin

    @test IdealGasThermo.saturation_vapor_pressure(273.16) ≈ 611.38 rtol = 1e-3
    @test IdealGasThermo.saturation_vapor_pressure(20 + 273.15) ≈ 2333.44 rtol = 1e-3
    @test IdealGasThermo.saturation_vapor_pressure(100 + 273.15) ≈ 104077 rtol = 1e-3

    @test IdealGasThermo.specific_humidity(1.0, 298.15, 101325.0) ≈ 0.0194 rtol = 1e-3
    @test IdealGasThermo.relative_humidity(0.0194, 298.15, 101325.0) ≈ 0.9996 rtol = 1e-3

    DA = IdealGasThermo.generate_humid_air(0.0)
    @test DA.composition == DryAir.composition

    RH = 0.1
    WetAir = IdealGasThermo.generate_humid_air(RH)
    SH = IdealGasThermo.specific_humidity(WetAir)

    @test SH ≈ IdealGasThermo.specific_humidity(RH, IdealGasThermo.Tstd, IdealGasThermo.Pstd)

end
