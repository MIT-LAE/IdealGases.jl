@testset "mix. comp." begin
    gas = Gas()
    Xair = Dict(
        "N2" => 0.78084,
        "Ar" => 0.009365,
        "Air" => 0.0,
        "H2O" => 0.0,
        "CO2" => 0.000319,
        "O2" => 0.209476,
    )

    Yair = Dict(
        "O2" => 0.231416,
        "Ar" => 0.012916,
        "Air" => 0.0,
        "H2O" => 0.0,
        "CO2" => 0.000484688,
        "N2" => 0.755184,
    )

    gas.X = IdealGases.Xair

    for (key, val) in gas.Ydict
        if val != 0
            @test Yair[key] ≈ val atol = 1e-6
        end
    end
    Air = species_in_spdict("Air")
    @test gas.MW ≈ Air.MW

end

@testset "mix. thermo" begin
    gas = Gas()
    gas.X = Xair = IdealGases.Xair

    Air = species_in_spdict("Air")

    # Low temp:
    gas.T = T = IdealGases.Tstd
    gas.P = P = IdealGases.Pstd

    @test gas.cp ≈ IdealGases.Cp(T, Air) rtol = 1e-7
    @test gas.h ≈ IdealGases.h(T, Air) rtol = 1e-7
    @test gas.s ≈ IdealGases.s(T, P, Air) rtol = 1e-7

    # High temp:
    gas.T = T = 20 * T
    gas.P = P = 20 * P
    @test gas.cp ≈ IdealGases.Cp(T, Air) rtol = 1e-7
    @test gas.h ≈ IdealGases.h(T, Air) rtol = 1e-7
    @test gas.s ≈ IdealGases.s(T, P, Air) rtol = 1e-7


    # Temp range test from https://cearun.grc.nasa.gov/cgi-bin/ThermoBuild/properties-3.pl
    #     THERMODYNAMIC FUNCTIONS CALCULATED FROM COEFFICIENTS FOR Air             

    #     T         Cp        H-H298        S      -(G-H298)/T      H        delta Hf     log K
    #   deg-K    J/mol-K      kJ/mol     J/mol-K     J/mol-K      kJ/mol      kJ/mol

    #    200      29.034      -2.852     187.221     201.481      -2.978      -0.125      0.2789
    #    300      29.105       0.054     199.002     198.823      -0.072      -0.126      0.2680
    #    400      29.355       2.975     207.404     199.967       2.849      -0.126      0.2625
    #    500      29.821       5.932     214.001     202.137       5.807      -0.126      0.2592
    #    600      30.442       8.944     219.491     204.584       8.819      -0.126      0.2570
    #    700      31.135      12.023     224.235     207.060      11.897      -0.126      0.2555
    #    800      31.825      15.171     228.438     209.474      15.046      -0.126      0.2543
    #    900      32.467      18.386     232.224     211.795      18.261      -0.126      0.2534
    #   1000      33.050      21.662     235.676     214.013      21.537      -0.126      0.2527
    #   1100      33.571      24.994     238.851     216.129      24.869      -0.126      0.2521
    #   1200      34.023      28.374     241.792     218.146      28.249      -0.126      0.2516
    #   1300      34.419      31.797     244.531     220.072      31.671      -0.127      0.2511
    #   1400      34.767      35.257     247.095     221.911      35.131      -0.127      0.2508
    #   1500      35.076      38.749     249.504     223.671      38.623      -0.127      0.2505
    #   1600      35.352      42.271     251.777     225.358      42.145      -0.127      0.2502
    #   1700      35.600      45.818     253.928     226.976      45.693      -0.127      0.2499
    #   1800      35.825      49.390     255.969     228.530      49.264      -0.127      0.2497
    #   1900      36.029      52.983     257.911     230.026      52.857      -0.127      0.2495
    #   2000      36.216      56.595     259.764     231.467      56.470      -0.127      0.2493

    Cp = [
        29.034,
        29.105,
        29.355,
        29.821,
        30.442,
        31.135,
        31.825,
        32.467,
        33.050,
        33.571,
        34.023,
        34.419,
        34.767,
        35.076,
        35.352,
        35.600,
        35.825,
        36.029,
        36.216,
    ]
    H = [
        -2.978,
        -0.072,
        2.849,
        5.807,
        8.819,
        11.897,
        15.046,
        18.261,
        21.537,
        24.869,
        28.249,
        31.671,
        35.131,
        38.623,
        42.145,
        45.693,
        49.264,
        52.857,
        56.470,
    ]
    S = [
        187.221,
        199.002,
        207.404,
        214.001,
        219.491,
        224.235,
        228.438,
        232.224,
        235.676,
        238.851,
        241.792,
        244.531,
        247.095,
        249.504,
        251.777,
        253.928,
        255.969,
        257.911,
        259.764,
    ]

    @testset "Temp range" begin
        Trange = 200.0:100.0:2000.0
        gas.P = IdealGases.Pstd
        Trange, cp_array, h_array, 𝜙_array, s_array = IdealGases.thermo_table(gas, Trange)
        for i in eachindex(Trange)
            @test cp_array[i] .* gas.MW ./ 1e3 ≈ Cp[i] atol = 1e-3
            @test h_array[i] .* gas.MW ./ 1e6 ≈ H[i] atol = 1e-3
            @test s_array[i] .* gas.MW ./ 1e3 ≈ S[i] atol = 1e-3
        end
    end

end
