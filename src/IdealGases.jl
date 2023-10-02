# """
# Thermally-perfect gas thermodynamics based on NASA polynomials
# """
module IdealGases

const __Gasroot__ = dirname(@__DIR__)
using LinearAlgebra
using StaticArrays
using Printf

export Gas, set_h!, set_hP!, set_TP!, set_Δh!

include("constants.jl")
include("species.jl")
export species

include("readThermo.jl")
include("Gas.jl")
include("combustion.jl")

include("io.jl")
include("utils.jl")
export X2Y, Y2X
include("idealgasthermo.jl")




end