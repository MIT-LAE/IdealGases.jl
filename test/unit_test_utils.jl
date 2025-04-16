@testset "utils" begin
    X = rand(0.0:0.001:1.0, IdealGasThermo.Nspecies)
    names = IdealGasThermo.spdict.name
    # Test dict/array conversions
    Xdict = Dict(zip(names, X))
    X1 = IdealGasThermo.Xidict2Array(Xdict)
    for i in eachindex(X)
        @test X1[i] == X[i]
    end
    IdealGasThermo.Xidict2Array!(Xdict, X1)
    X = X / sum(X)
    for i in eachindex(X)
        @test X1[i] == X[i]
    end
    # Test mass/molar fraction conversions
    Y = IdealGasThermo.X2Y(X)
    X2 = IdealGasThermo.Y2X(Y)
    for i in eachindex(X)
        @test X[i] ≈ X2[i]
    end
end
