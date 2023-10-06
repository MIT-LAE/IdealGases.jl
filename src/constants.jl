const Runiv = 8.3145 # J/K/mol
const Pstd = 101325.0 # Pa
const Tstd = 298.15 # K
const ϵ = sqrt(eps()) #standard tolerance 

# Air composition
const Xair = Dict(  
"N2"  => 0.78084,
"Ar"  => 0.009365,
"Air" => 0.0,
"H2O" => 0.0,
"CO2" => 0.000319,
"O2"  => 0.209476)