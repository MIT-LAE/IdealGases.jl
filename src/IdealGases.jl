# """
# Thermally-perfect gas thermodynamics based on NASA polynomials
# """
module IdealGases

const __Gasroot__ = dirname(@__DIR__)
const default_thermo_path = joinpath(__Gasroot__, "data/thermo.inp")
using LinearAlgebra
using StaticArrays
using Printf

export Gas, set_h!, set_hP!, set_TP!, set_Δh!

include("constants.jl")
include("species.jl")
export AbstractSpecies, species, composite_species, generate_composite_species

include("readThermo.jl")
export readThermo, species_in_spdict
include("Gas.jl")
include("combustion.jl")

include("io.jl")
export print_thermo_table
include("utils.jl")
export X2Y, Y2X
include("idealgasthermo.jl")




end