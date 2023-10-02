@testset "read thermo" begin

species_dict = readThermo(IdealGases.default_thermo_path)
Air = species_dict[findfirst(x->x=="Air", species_dict.name)]

@test Air.MW == 28.9651159
@test Air.Hf == -125.530
@test Air.alow[1] == 1.009950160e+04
@test Air.Tmid == 1000.0
@test Air.ahigh[9] == -8.147411905


end