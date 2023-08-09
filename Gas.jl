# """
# Thermally-perfect gas thermodynamics based on NASA polynomials
# """
# module IdealGas

# using NLsolve
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf

include("constants.jl")
include("readThermo.jl")

"""
    Gas

A type that represents an ideal gas that is calorically perfect 
"""
mutable struct Gas
   P::Float64 # [Pa]
   T::Float64 # [K]
   Tarray::MVector{8, Float64} # Temperature array to make calcs allocation free

   cp::Float64 #[J/mol/K]
   h::Float64  #[J/mol]
   s::Float64  #[J/mol/K]
   Y::MVector{length(spdict), Float64} # Mass fraction of species
   MW::Float64 # Molecular weight [g/mol]
end

"""
    Gas(Y)

Constructs `Gas` with given composition `Y`

"""
function Gas(Y)
   Gas(Pstd, Tstd, Tarray(Tstd), 0.0, 0.0, 0.0, Y, 28.965)
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
    h(Tstd, Air),
    s(Tstd, Pstd, Air),
   [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], Air.MW)

end

# Overload Base.getproperty for convinence
function Base.getproperty(gas::Gas, s::Symbol)
   if s === :h_T
      return getfield(gas, :cp)
   elseif s === :s_T
      return getfield(gas, :cp)/getfield(gas, :T)
   elseif s === :hs
      return [getfield(gas, :h), getfield(gas, :s)]
   elseif s === :TP
      return [getfield(gas, :T), getfield(gas, :P)]
   else
      return getfield(gas, s)
   end
end
"""
   show(io::IO, gas::Gas)

Pretty print for Real gases
"""
function Base.show(io::IO, gas::Gas)
   @printf(io, "Ideal Gas at\n%3s = %8.3f K\n%3s = %8.3f kPa\n%3s = %8.3f J/K/mol\n%3s = %8.3f kJ/mol\n%3s = %8.3f kJ/K/mol",
     "T", gas.T, "P", gas.P/1000.0, "cp", gas.cp, "h", gas.h/1000.0, "s", gas.s/1000.0)
   println("\n\nwith composition:")
   composition(gas,io)
end

"""
    composition(gas::Gas, io::IO=stdout)

Prints out the composition (Y_i) of the gas
"""
function composition(gas::Gas, io::IO=stdout)
   divider = "-"^(8*2+4+9)
   @printf(io, "%s\n", divider)
   @printf(io, "%8s  %8s  %9s\n", "Species", "Yᵢ", "MW[g/mol]")
   @printf(io, "%s\n", divider)
   for (name, Yi, mw) in zip(spdict.name, gas.Y, spdict.MW)
      if Yi != 0
         @printf(io, "%8s  %8.3f  %9.3f\n", name, Yi, mw)
      end
   end
   @printf(io, "%s\n", divider)
   @printf(io, "%8s  %8.3f  %9.3f\n", "ΣYᵢ", sum(gas.Y), gas.MW)
end


"""
    Base.setproperty!(gas::Gas, s::Symbol, val)

"""
function Base.setproperty!(gas::Gas, s::Symbol, val)
   if s === :T
      setfield!(gas, :T, val) # first set T
      setfield!(gas, :Tarray, Tarray!(val, getfield(gas, :Tarray))) # update Tarray
      TT = getfield(gas, :Tarray) # Just convinence
      # Next set the cp, h and s of the gas
      ## Get the right coefficients (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
      if val<1000.0
         A = view(spdict.alow, :)
      else
         A = view(spdict.ahigh, :)
      end   
      ## Initialize temporary vars
      cptemp = 0.0
      htemp  = 0.0
      stemp  = 0.0
      
      P = getfield(gas, :P)
      Y = getfield(gas, :Y)
      # Go through every species where mass fraction is not zero
      @views for (Yᵢ,a) in zip(Y, A)
         if Yᵢ != 0.0
            cptemp = cptemp + Yᵢ * Cp(TT, a)
             htemp = htemp  + Yᵢ * h(TT, a)
             stemp = stemp  + Yᵢ * (𝜙(TT, a) - Runiv*log(P/Pstd))
         end
      end
   
      setfield!(gas, :cp, cptemp)
      setfield!(gas, :h, htemp)
      setfield!(gas, :s, stemp)

   elseif s === :P
      setfield!(gas, :P, val)
      TT = view(getfield(gas, :Tarray), :) # Just convinence
      # Next set s of the gas
      ## Get the right coefficients (assumes Tmid is always 1000.0. Check performed in readThermo.jl.):
      if val<1000.0
         A = view(spdict.alow, :)
      else
         A = view(spdict.ahigh, :)
      end   
      ## Initialize temporary vars
      stemp  = 0.0
      
      P = val
      Y = view(getfield(gas, :Y), :)
      # Go through every species where mass fraction is not zero
      @views for (Yᵢ,a) in zip(Y, A)
         if Yᵢ != 0.0
            stemp = stemp  + Yᵢ * (𝜙(TT, a) - Runiv*log(P/Pstd))
         end
      end

      setfield!(gas, :s, stemp)

   elseif s === :Y # directly set mass fractions Y
      if typeof(val) === Array{Float64, 1}
         # If array directly store in Y
         setfield!(gas, :Y, val) 
      elseif typeof(val) <: Dict
         # If dict provided set each species in the right order
         names = spdict.name
         Y = zeros(MVector{length(names)})
         for (key,value) in val
            index = findfirst(x->x==key, names)
            Y[index] = value
         end
         setfield!(gas, :Y, Y)
      end
      # Update the MW of the gas mixture
      setfield!(gas, :MW, MW(gas))
   end
   # Note: intentionally not including other variables to prevent users from trying to directly set h, s, cp, MW etc.

end

"""
Function to create the required temperature array
"""
function Tarray(T)
   return [T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
end

# function Tarray2(T)
#    TT = zeros(Float64, 8)
#    TT[1:7] = [T^i for i in range(-2, stop=4)]
#    TT[8] = log(T)
#    return TT
# end

"""
In place Tarray update that returns
[T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
"""
function Tarray!(T, TT)
   TT[1] = T^-2    #T^-2
   TT[2] = TT[1]*T #T^-1
   TT[3] = 1.0     #T^0
   TT[4] = T       #T^1
   TT[5] = T*T     #T^2
   TT[6] = T*TT[5] #T^3
   TT[7] = T*TT[6] #T^4
   TT[8] = log(float(T))
   return TT
end


"""
Calculates cp of the given species in J/K/mol
(This is a completely non-allocating operation.)
"""
@views function Cp(Tarray, a)
   #  Cp_R = dot(view(a, 1:7), view(Tarray, 1:7))
    Cp_R = dot(a[1:7], Tarray[1:7])
    Cp = Cp_R*Runiv
    return Cp #J/K/mol
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
Calculates cp of a mixture specified by the mass fraction in `gas`
"""
@views function Cp(T, g::Gas)
   g.T = T
   return g.cp
end

function Cp(g::Gas)
   Cp(g.T, g)
end

"""
Calculates mean molecular weight
"""
@views function MW(g::Gas)
   MW = dot(g.Y, spdict.MW)
   return MW
end

"""
Calculates h of the given **species** in J/mol
Calcualted by:
H0/RT = -a1*T^-2 + a2*T^-1*ln(T) + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
      = -a1*T₁   + a2*T₂*T₈      + a3 + a4*T₄/2 + a5*T₅/3  + a6*T₆/4  + a7*T₇/5  + a₈*T₂
"""
function h(TT, a)
    h_RT  = -a[1]*TT[1] + 
             a[2]*TT[8]*TT[2] + 
             a[3] + 
         0.5*a[4]*TT[4] + 
             a[5]*TT[5]/3.0 + 
        0.25*a[6]*TT[6] + 
        0.20*a[7]*TT[7] + 
             a[8]*TT[2]

    h = h_RT*TT[4]*Runiv
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
Calculates h of a given **mixture** in J/mol where species mass fractions \\math{Y_i} 
is calculated from the supplied Gas instance
"""
function h(T, g::Gas)
   H = 0.0
   g.T = T
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   
   for (key,Yᵢ) in g.Y
      a = getfield(spd[key], s)
      H = H + Yᵢ * h(g.Tarray, a)
   end
   return H
end
function h(g::Gas)
   h(g.T,g)
end

"""
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
Calculates the entropy complement function 𝜙=∫(cₚ/T)dT of the given **mixture** in J/K/mol
This is calculated at standard state. Tref = 298.15 K, Pref = 101325 Pa.
"""
function 𝜙(T, g::Gas)
   S = 0.0
   g.T = T
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end

   for (key,Yᵢ) in g.Y
      a = getfield(spd[key], s)
      S = S + Yᵢ * 𝜙(g.Tarray, a)
   end
   return S
end
function 𝜙(g::Gas)
   𝜙(g.T, g)
end


"""
Returns standard state sᵒ based on the reference point defined at
Tref = 298.15 K
Pref = 101325 Pa

using the entropy complement function and
the entropy change due to pressure.
Δs can then be defined as sᵒ - sᵒ(Tref, Pref) = sᴼ - 𝜙
"""
function s(T, P, gas::Gas)
   Pref = 101325.0 
   gas.T = T
   sᵒ =  𝜙(gas) - Runiv*log(P/Pref)
   return sᵒ
end

"""
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
"""
function set_h!(gas::Gas, hspec::Float64)
   T = gas.T
   dT = T
   while abs(dT) > ϵ
      h = gas.h
      res = h - hspec # Residual
      res_t = gas.cp  # ∂R/∂T = ∂h/∂T = cp
      dT = -res/res_t # Newton step
      T = T + dT
      gas.T = T
   end
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


# # Specific functions for gas Compression
# PR = 10
# p2, T2 = 101325, 298.15
# ηp = 0.90

# """ 
# Adiabatic compression given the 
# compression pressure ratio (`PR`), the initial pressure (`p`)
# and initial temperature (`T`).

# Returns `Tfinal` and `pfinal`
# """
# function compress(PR, p, T)
#    Tfinal = T * PR^(ℜ/cp(T,Air))

#    for i in 1:10
#       Res  = (𝜙(Tfinal, Air) - 𝜙(T, Air))/ℜ - log(PR)
#       Res′ = cp(Tfinal,Air)/ℜ/Tfinal
#       dT  = Res/Res′
#       Tfinal = Tfinal - dT
#       # println(Tfinal)
#       if abs(dT) < ϵ
#          break
#       end
#    end

#    return Tfinal, p*PR

# end
# """
# Adiabatic with NL solve
# i.e. find x such that F(x)=0
# """
# T = 298.15
# p = 101325.
# PR = 2.0
# function f(x)
#    s(T,p,Air) - s(x[1],p*PR,Air)
# end

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
   
   while abs(dT)>ϵ
      ## Original approach by M. Drela using entropy complement
      # res  = (𝜙(Tfinal, Air) - s)/Runiv - log(PR)/ηp
      # res_dT = cp(Tfinal,Air)/Runiv/Tfinal
      ## Modified approach using pressure dependent entropy
      res  = (gas.s - s0)/Runiv + (log(PR) - log(PR)/ηp)
      res_dT = gas.s_T/Runiv
      dT  = - res/res_dT

      Tfinal = Tfinal + dT
      gas.T = Tfinal
      # println("$i: $Tfinal $dT")
   end

   return gas

end

"""
    thermo_table(gas::Gas, Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0)

Quickly generate a table of cp, h and s for a gas
"""
function thermo_table(gas::Gas, Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0)
   Trange = range(Tstart, Tend, step=Tinterval)
   composition(gas)
   println(" ")
   divider = "-"^(4+8+12*3+4)
   @printf("%4s %8s %12s %12s %12s\n",
          "i",  "T[K]", "cp[J/K/mol]", "h[kJ/mol]", "s[kJ/K/mol]")
   println(divider)
   for (i,T) in enumerate(Trange)
      gas.T = T
      @printf("%4d %8.2f %12.4f %12.4f %12.4f\n",
      i,  gas.T, gas.cp, gas.h/1000.0, gas.s/1000.0)
   end
end

Y = Dict(  
"N2"  => 0.78084,
"Ar"  => 0.009365,
"Air" => 0.0,
"H2O" => 0.0,
"CO2" => 0.000319,
"O2"  => 0.209476)

gas = Gas()

