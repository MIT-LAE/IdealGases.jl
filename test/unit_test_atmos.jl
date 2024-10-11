@testset "atmosphere" begin
    #Data from US Standard Atmosphere tables for set geometric altitude
    Zs = [8, 18, 28, 38, 48, 58, 68, 78] * 1e3
    Ts = [236.215, 216.65, 224.527, 244.818, 270.65, 252.518, 225.065, 202.541]

    Ps = [35651, 7565.2, 1616.1, 377.13, 102.29, 28.723, 7.0529, 1.4673]

    for (i, Z) in enumerate(Zs)
        gas = IdealGases.standard_atmosphere(Z)

        @test gas.T ≈ Ts[i] rtol = 1e-4
        @test gas.P ≈ Ps[i] rtol = 1e-4
    end

end
