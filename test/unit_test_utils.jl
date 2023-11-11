@testset "utils" begin
    X = rand(0.0:0.001:1.0, IdealGases.Nspecies)
    names = IdealGases.spdict.name
    # Test dict/array conversions
    Xdict = Dict(zip(names, X))
    X1 = IdealGases.Xidict2Array(Xdict)
    for i in eachindex(X)
        @test X1[i] == X[i]
    end
    IdealGases.Xidict2Array!(Xdict, X1)
    X = X/sum(X)
    for i in eachindex(X)
        @test X1[i] == X[i]
    end
    # Test mass/molar fraction conversions
    Y = IdealGases.X2Y(X)
    X2 = IdealGases.Y2X(Y)
    for i in eachindex(X)
        @test X[i] ≈ X2[i]
    end
end