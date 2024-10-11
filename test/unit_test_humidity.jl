@testset "humidity" begin

    @test IdealGases.saturation_vapor_pressure(273.16) ≈ 611.38 rtol = 1e-3
    @test IdealGases.saturation_vapor_pressure(20 + 273.15) ≈ 2333.44 rtol = 1e-3
    @test IdealGases.saturation_vapor_pressure(100 + 273.15) ≈ 104077 rtol = 1e-3

    @test IdealGases.specific_humidity(1.0, 298.15, 101325.0) ≈ 0.0194 rtol = 1e-3
    @test IdealGases.relative_humidity(0.0194, 298.15, 101325.0) ≈ 0.9996 rtol = 1e-3

    DA = IdealGases.generate_humid_air(0.0)
    @test DA.composition == DryAir.composition

    RH = 0.1
    WetAir = IdealGases.generate_humid_air(RH)
    SH = IdealGases.specific_humidity(WetAir)

    @test SH ≈ IdealGases.specific_humidity(RH, IdealGases.Tstd, IdealGases.Pstd)

end
