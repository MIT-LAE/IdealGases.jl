using IdealGases

Xair = Dict(  
"N2"  => 0.78084,
"Ar"  => 0.009365,
"H2O" => 0.0,
"CO2" => 0.000319,
"O2"  => 0.209476)

Xox = IdealGases.Xidict2Array(Dict([("O2", 1.0)])) #Mole fraction
Yox = X2Y(Xox) 

Xf = IdealGases.Xidict2Array(Dict([("H2", 1.0)])) #Mole fraction
Yf = X2Y(Xf) 

fuel_sps = species_in_spdict("H2")
gas_sps = generate_composite_species(Xox)

prod_species = generate_composite_species(IdealGases.Xidict2Array(Dict([("H2O", 1.0)])))

gas_ox = Gas(Yox) 
gasf = Gas(Yf) 

FAR = 0.1259972248959336
gas2 = IdealGases.fuel_combustion(gas_ox, "H2", 298.15, 0.1259972248959336)

Xw = IdealGases.Xidict2Array(Dict([("H2O", 1.0)])) #Mole fraction
Yw = X2Y(Xw) 
gas_w = Gas(Yw) 

if "Air" in keys(gas_ox.Xdict)
    Xin = Xair
else
    Xin = gas_ox.Xdict
end
gas_sps = generate_composite_species(IdealGases.Xidict2Array(Xin))

Xprod_dict =  IdealGases.vitiated_mixture(fuel_sps, gas_sps, 0.1259972248959336, 1.0)

#Initialize output 
gas_prod = Gas1D(prod_species)