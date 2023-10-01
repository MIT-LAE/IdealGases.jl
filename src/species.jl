"""
species is a structure that holds the NASA 9 polynomial coefficients `alow` and `ahigh` 
for the two temprature regions separated by `Tmid` 
(here we only work with temperature less than 6000 K so typically only 2 T intervals required)
the molecular weight `MW` and the heat of formation `Hf` (J/mol) for a given chemical species (at 298.15 K).

See https://shepherd.caltech.edu/EDL/PublicResources/sdt/formats/nasa.html for typical data format
"""
struct species
    name::String
    Tmid::Float64
    alow::Array{Float64, 1}
    ahigh::Array{Float64, 1}
    MW::Float64
    Hf::Float64
end
