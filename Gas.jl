# """
# Thermally-perfect gas thermodynamics based on NASA polynomials
# """
# module IdealGas

# using NLsolve
using LinearAlgebra
using BenchmarkTools, ProfileView
using StaticArrays
using Printf

include("constants.jl")
include("readThermo.jl")
include("combustion.jl")

"""
    Gas

A type that represents an ideal gas that is calorically perfect 
"""
mutable struct Gas
   P::Float64 # [Pa]
   T::Float64 # [K]
   Tarray::MVector{8, Float64} # Temperature array to make calcs allocation free

   cp::Float64 #[J/mol/K]
   cp_T::Float64 # derivative dcp/dT
   h::Float64  #[J/mol]
   ϕ::Float64  #[J/mol/K] Entropy complement fn ϕ(T) = ∫ cp(T)/T dT from Tref to T
   
   Y::MVector{length(spdict), Float64} # Mass fraction of species
   MW::Float64 # Molecular weight [g/mol]
end

"""
    Gas(Y)

Constructs `Gas` with given composition `Y`

"""
function Gas(Y::AbstractVector)
   gas = Gas(); gas.Y = convert(Vector{Float64}, Y)
   return gas
end

"""
    Gas()

Constructor that returns a `Gas` type representing 
Air at standard conditions

See also [`Gas`](@ref).

# Examples
```julia-repl
julia> Gas()
Ideal Gas at
  T =  298.150 K
  P =  101.325 kPa
 cp =   29.102 J/K/mol
  h =   -0.126 kJ/mol
  s =    0.199 kJ/K/mol

with composition:
-----------------------------
 Species        Yᵢ  MW[g/mol]
-----------------------------
     Air     1.000     28.965
-----------------------------
     ΣYᵢ     1.000     28.965
```
"""
function Gas()
   Air = spdict[findfirst(x->x=="Air", spdict.name)]

   Gas(Pstd, Tstd, Tarray(Tstd),
    Cp(Tstd, Air), 
    (Cp(Tstd + 1.0, Air) - Cp(Tstd - 1.0, Air) /2.0), #finite diff dCp/dT
    h(Tstd, Air),
    s(Tstd, Pstd, Air),
   [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], Air.MW)

end

# Overload Base.getproperty for convinence
function Base.getproperty(gas::Gas, sym::Symbol)
   if sym === :h_T # dh/dT
      return getfield(gas, :cp)
   elseif sym === :s_T # ∂h/∂T
      return getfield(gas, :cp)/getfield(gas, :T)
   elseif sym === :hs
      return [getfield(gas, :h), getfield(gas, :s)]
   elseif sym === :TP
      return [getfield(gas, :T), getfield(gas, :P)]
   elseif sym === :s
      return getfield(gas, :ϕ) - Runiv*log(getfield(gas,:P)/Pstd)
   else
      return getfield(gas, sym)
   end
end


"""
    Base.setproperty!(gas::Gas, s::Symbol, val)

"""
function Base.setproperty!(gas::Gas, sym::Symbol, val::Float64)
   ## Setting Temperature
   if sym === :T
      setfield!(gas, :T, val) # first set T
      setfield!(gas, :Tarray, Tarray!(val, getfield(gas, :Tarray))) # update Tarray
      TT = view(getfield(gas, :Tarray), :) # Just convinence
      # Next set the cp, h and s of the gas
      ## Get the right coefficients 
      ## (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
      if val<1000.0 #Value given is the desired temperature
         A = view(spdict.alow, :)
      else
         A = view(spdict.ahigh, :)
      end   
      ## Initialize temporary vars
      cptemp = 0.0
      htemp  = 0.0
      ϕtemp  = 0.0
      cp_Ttemp = 0.0

      P = getfield(gas, :P)
      Y = getfield(gas, :Y)
      # Go through every species where mass fraction is not zero
      @inbounds for (Yᵢ,a) in zip(Y, A)
         if Yᵢ != 0.0
            cptemp = cptemp + Yᵢ * Cp(TT, a)
             htemp = htemp  + Yᵢ * h(TT, a)
             ϕtemp = ϕtemp  + Yᵢ * 𝜙(TT, a)
             cp_Ttemp = cp_Ttemp + Yᵢ * dCpdT(TT, a) 
         end
      end
   
      setfield!(gas, :cp, cptemp)
      setfield!(gas, :h, htemp)
      setfield!(gas, :ϕ, ϕtemp)

   ## Setting Pressure
   elseif sym === :P
      setfield!(gas, :P, val)
      TT = view(getfield(gas, :Tarray), :) # Just convinence
      # Next set s of the gas
      ## Get the right coefficients 
      ## (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
      if TT[4]<1000.0 #TT[4] is T
         A = view(spdict.alow, :)
      else
         A = view(spdict.ahigh, :)
      end   
      ## Initialize temporary vars
      ϕtemp  = 0.0
      
      P = val
      Y = view(getfield(gas, :Y), :)
      # Go through every species where mass fraction is not zero
      @inbounds for (Yᵢ,a) in zip(Y, A)
         if Yᵢ != 0.0
            ϕtemp = ϕtemp  + Yᵢ * 𝜙(TT, a)
         end
      end

      setfield!(gas, :ϕ, ϕtemp)

   elseif sym ===:h
      set_h!(gas, val)
   elseif sym === :TP
      set_TP!(gas, val[1], val[2])
   else
      error("You tried setting gas.", sym, " to a ", typeof(val))
   end
   # Note: intentionally not including other variables to prevent 
   # users from trying to directly set h, s, cp, MW etc.
   nothing
end

function Base.setproperty!(gas::Gas, sym::Symbol, val::AbstractVector{Float64})
   if sym === :Y # directly set mass fractions Y
      n = length(getfield(gas, :Y))
      setfield!(gas, :Y, MVector{n}(val)) 
      setfield!(gas, :MW, MW(gas)) # Update the MW of the gas mixture

   else
      error("Only mass factions Y can be set with an array.",
      "\n       You tried to set gas.",sym," to a ",typeof(val))
   end
   nothing
end

function Base.setproperty!(gas::Gas, sym::Symbol, val::AbstractDict{String, Float64})
   if sym === :Y # directly set mass fractions Y
      # If dict provided set each species in the right order
      names = spdict.name
      Y = zeros(MVector{length(names), Float64})
      for (key,value) in val
         index = findfirst(x->x==key, names)
         Y[index] = value
      end
      setfield!(gas, :Y, Y)
      setfield!(gas, :MW, MW(gas)) # Update the MW of the gas mixture
   else
      error("Only mass factions Y can be set with a dict.",
      "\n       You tried to set gas.",sym," to a ",typeof(val))
   end
   nothing
end

include("io.jl")
include("utils.jl")

"""
Calculates cp of the given species in J/K/mol
(This is a completely non-allocating operation.)

Cp0/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
"""
@views function Cp(Tarray::AbstractVector{T}, a::AbstractVector{T}) where T
   #  Cp_R = dot(view(a, 1:7), view(Tarray, 1:7))
    Cp_R = dot(a[1:7], Tarray[1:7])
    Cp = Cp_R*Runiv
    return Cp #J/K/mol
end

"""
    dCpdT(TT::AbstractVector{T}, a::AbstractVector{T}) where T

Returns the derivative dcp/dT
dCp0/dT = R×(-2a1*T^-3 -a2*T^-2 + a4 + 2a5*T + 3a6*T^2 + 4a7*T^3)
"""
function dCpdT(TT::AbstractVector{T}, a::AbstractVector{T}) where T
   dcp_RdT = -2*a[1]*TT[1]*TT[2] -
             a[2]*TT[1] +
             a[4] +
           2*a[5]*TT[4] +
           3*a[6]*TT[5] +
           4*a[7]*TT[6]
   return dcp_RdT * Runiv
end
"""
Calculates cp for a **species** type in J/K/mol.
"""
function Cp(T, sp::species)
   TT = Tarray(T)
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   a = getfield(sp, s)
   Cp(TT, a)
end

"""
Calculates mean molecular weight
"""
@views function MW(g::Gas)
   MW = dot(g.Y, spdict.MW)
   return MW
end

"""
    h(TT::AbstractVector{type}, a::AbstractVector{type}) where type

Calculates h of the given **species** in J/mol
Calcualted by:
H0/RT = -a1*T^-2 + a2*T^-1*ln(T) + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
      = -a1*T₁   + a2*T₂*T₈      + a3 + a4*T₄/2 + a5*T₅/3  + a6*T₆/4  + a7*T₇/5  + a₈*T₂
"""
function h(TT::AbstractVector{type}, a::AbstractVector{type}) where type
    h_RT  = -a[1]*TT[1] + 
             a[2]*TT[8]*TT[2] + 
             a[3] + 
         0.5*a[4]*TT[4] + 
             a[5]*TT[5]/3.0 + 
        0.25*a[6]*TT[6] + 
        0.20*a[7]*TT[7] + 
             a[8]*TT[2]

    h = h_RT*TT[4]*Runiv # because TT[4] == T
    return h #J/mol
end

"""
Calculates h for a species
"""
function h(T, sp::species)
   TT = Tarray(T)
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   a = getfield(sp, s)
   h(TT, a)
end

"""
    𝜙(TT,a)

Calculates the entropy complement function 𝜙=∫(cₚ/T)dT in J/K/mol
This is calculated at standard state. Tref = 298.15 K, Pref = 101325 Pa.
```math
S0/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3.0 + a7*T^4/4 + b2 
     = -a1*T₁/2   - a2*T₂   + a3*T₈    + a4*T₄+ a5*T₅/2  + a6*T₆/3.0  + a7*T₇/4  + a₉   
```
"""
function 𝜙(TT,a)
    so_R = -0.5*a[1] * TT[1] - 
                a[2] * TT[2] + 
                a[3] * TT[8] + 
                a[4] * TT[4] + 
            0.5*a[5] * TT[5] +
                a[6] * TT[6]/3.0 + 
           0.25*a[7] * TT[7] + 
                a[9]

    so = so_R*Runiv
    return so #J/K/mol
end

"""
    𝜙(g::Gas)

Calculates the entropy complement function 𝜙=∫(cₚ/T)dT of the 
given **mixture** in J/K/mol
This is calculated at standard state. Tref = 298.15 K, Pref = 101325 Pa.
"""
function 𝜙(g::Gas)
   ϕ = 0.0
   if gas.T<1000.0
      A = view(spdict.alow, :)
   else
      A = view(spdict.ahigh, :)
   end

   for (Yᵢ, a) in zip(g.Y, A)
      ϕ = ϕ + Yᵢ * 𝜙(g.Tarray, a)
   end
   return ϕ
end

"""
    s(T, P, sp::species)

Calculates s for a species
"""
function s(T, P, sp::species)
   TT = Tarray(T)
   Pref = 101325
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   a = getfield(sp, s)
   sᵒ = 𝜙(TT, a) - Runiv*log(P/Pref)
   return sᵒ
end

"""
    set_h!(gas::Gas, hspec::Float64)

Calculates gas temperature for a specified enthalpy

# Examples
```julia-repl
julia> gas = Gas();
julia> set_h!(gas, 0.0)
Ideal Gas at
  T =  302.463 K
  P =  101.325 kPa
 cp =   29.108 J/K/mol
  h =    0.000 kJ/mol
  s =    0.199 kJ/K/mol

with composition:
-----------------------------
 Species        Yᵢ  MW[g/mol]
-----------------------------
     Air     1.000     28.965
-----------------------------
     ΣYᵢ     1.000     28.965
```
"""
function set_h!(gas::Gas, hspec::Float64)
   
   T = gas.T
   dT = T
   for i in 1:20 # abs(dT) > ϵ
      res = gas.h - hspec # Residual
      res_t = gas.cp  # ∂R/∂T = ∂h/∂T = cp
      dT = -res/res_t # Newton step
      
      if abs(dT) ≤ ϵ
         break
      end

      T = T + dT
      gas.T = T
   end

   if abs(dT) > ϵ
      error("Error: `set_h!` did not converge:\ngas=", print(gas),
       "\n\nabs(dT) = ", abs(dT), " > ϵ (", ϵ, ")")
   end

   return gas
end
"""
    set_Δh!(gas::Gas, Δhspec::Float64)

Sets the gas temperature based on a specified change in enthalpy (Δh) [J/mol]
"""
function set_Δh!(gas::Gas, Δhspec::Float64)
   hf = gas.h + Δh
   set_h!(gas, hf)
   return gas
end
"""
    set_hP!(gas::Gas, hspec::Float64, P::Float64)

Calculates state of the gas given enthalpy and pressure (h,P)
"""
function set_hP!(gas::Gas, hspec::Float64, P::Float64)
   set_h!(gas,hspec)
   gas.P = P
   return gas
end
"""
    set_TP!(gas::Gas, T::Float64, P::Float64)

Calculates state of the gas given Temperature and pressure (T,P)
in K and Pa respectively.

# Examples
```julia-repl
julia> gas = Gas(); # Create an ideal gas consisting of air at std. conditions
julia> set_TP!(gas, 298.15*2, 101325.0*2)
Ideal Gas at
  T =  596.300 K
  P =  202.650 kPa
 cp =   30.418 J/K/mol
  h =    8.706 kJ/mol
  s =    0.214 kJ/K/mol

with composition:
-----------------------------
 Species        Yᵢ  MW[g/mol]
-----------------------------
     Air     1.000     28.965
-----------------------------
     ΣYᵢ     1.000     28.965
```

"""
function set_TP!(gas::Gas, T::Float64, P::Float64)
   gas.T = T
   gas.P = P
   return gas
end


"""
    compress(gas::Gas, PR::Float64, ηp::Float64=1.0,)

Compression with polytropic efficiency
"""
function compress(gas::Gas, PR::Float64, ηp::Float64=1.0,)

   T0 = gas.T
   s0 = gas.s
   P0 = gas.P

   Tfinal = T0 * PR^(Runiv/gas.cp/ηp)
   Pfinal = P0*PR
   dT = Tfinal
   gas.P = Pfinal
   gas.T = Tfinal
   
   for i in 1:20# abs(dT)>ϵ
      ## Original approach by M. Drela using entropy complement
      # res  = (𝜙(Tfinal, Air) - s)/Runiv - log(PR)/ηp
      # res_dT = cp(Tfinal,Air)/Runiv/Tfinal
      ## Modified approach using pressure dependent entropy
      res  = (gas.s - s0)/Runiv + (log(PR) - log(PR)/ηp)
      res_dT = gas.s_T/Runiv
      dT  = - res/res_dT

      if abs(dT) ≤ ϵ
         break
      end

      Tfinal = Tfinal + dT
      gas.T = Tfinal
      # println("$i: $Tfinal $dT")
   end

   if abs(dT) > ϵ
      error("Error: Compress did not converge:\ngas=", print(gas),
       "\n\nabs(dT) = ", abs(dT), " > ϵ (", ϵ, ")")
   end

   return gas

end

Y = Dict(  
"N2"  => 0.78084,
"Ar"  => 0.009365,
"Air" => 0.0,
"H2O" => 0.0,
"CO2" => 0.000319,
"O2"  => 0.209476)

gas = Gas()

