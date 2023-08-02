using NLsolve
using LinearAlgebra
using BenchmarkTools

const Runiv = 8.3145 # J/K/mol

Y = Dict()
Y = Dict(  
"N2"  => 0.78084,
"Ar"  => 0.009365,
"Air" => 0.0,
"H2O" => 0.0,
"CO2" => 0.000319,
"O2"  => 0.209476)


mutable struct Gas
   T::Float64
   Tarray::Array{Float64, 1} # Temperature array of gas
   Y # Mass fraction of species
end
# Convinence constructors:
function Gas(T, Y)
   Gas(T, Tarray(T), Y)
end
function Gas(Y)
   Gas(298.0, Y)
end
function Gas()
   Gas(298.0, Dict("Air"=>1.0))
end

# Automatically calculates the Tarray if T is set
function Base.setproperty!(gas::Gas, s::Symbol, val)
   if s === :T
      setfield!(gas, :T, val) # first set T
      setfield!(gas, :Tarray, Tarray!(val, getfield(gas, :Tarray))) # update Tarray
   elseif s === :Y # directly set mass fractions Y
      setfield!(gas, :Y, val) 
   elseif s === :Tarray
      setfield!(gas,:Tarray, val)
      setfield!(gas, :T, val[4])
   end

end

#Read species data from thermo.inp
spd = readThermo("thermo.inp")

Air = gas(Y)


ϵ = 1e-12

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
function Tarray!(T, TT)
   TT[1] = T^-2    #T^-2
   TT[2] = TT[1]*T #T^-1
   TT[3] = 1.0     #T^1
   TT[4] = T       #T^0
   TT[5] = T*T     #T^2
   TT[6] = T*TT[5] #T^3
   TT[7] = T*TT[6] #T^4
   TT[8] = log(float(T))
   return TT
end


"""
Calculates cp of the given species in J/K/mol
"""
function cp(T, sp::species)
    Tarray = [T^i for i in range(-2, stop = 4)]
    if T<1000.0
      a = sp.a_dict[200.0]
    else
      a = sp.a_dict[1000.0]
    end
    Cp_R   = a[1:7]' * Tarray
    Cp = Cp_R*ℜ
    return Cp #J/K/mol
end
"""
Calculates cp of a mixture specified by the mass fraction in `gas`
"""
function cp(T, g::gas)
   Cp = 0
   for (key,val) in g.Y
      Cp = Cp + val * cp(T, spd[key])
   end
   return Cp
end

"""
Calculates mean molecular weight
"""
function MW(g::gas)
   MW = 0
   for (key,val) in g.Y
      MW = MW + val*spd[key].MW
   end
   return MW/1000
end
"""
Calculates mean molecular weight
"""
function MW(g::species)
   MW = g.MW
   return MW/1000
end

"""
Calculates h of the given **species** in J/mol
"""
function h(T, sp::species)
   if T<1000.0
      a = sp.a_dict[200.0]
   else
      a = sp.a_dict[1000.0]
   end
    h_R_T = -a[1] * T^-2 + a[2]*log(T) * T^-1 + a[3] * T^0 + a[4]/2 * T + a[5]/3 * T^2 + a[6]/4 * T^3 + a[7]/5 * T^4 + a[8]/T
    h = h_R_T*T*ℜ
    return h #J/mol
end
"""
Calculates h of a given **mixture** in J/mol
"""
function h(T, g::gas)
   H = 0
   for (key,val) in g.Y
      H = H + val * h(T, spd[key])
   end
   return H
end

"""
Calculates the entropy complement function 𝜙=∫(cₚ/T)dT of the given **species** in J/K/mol
"""
function 𝜙(T::Float64,sp::species)
   if T<1000.0
      a = sp.a_dict[200.0]
   else
      a = sp.a_dict[1000.0]
   end

    so_R = -a[1]/2 * T^-2 - a[2]* T^-1 + a[3] * log(T) + a[4] * T + a[5]/2 * T^2 + a[6]/3 * T^3 + a[7]/4 * T^4 + a[9]
    so = so_R*ℜ
    return so #J/K/mol
end
"""
Calculates the entropy complement function 𝜙=∫(cₚ/T)dT of the given **mixture** in J/K/mol
"""
function 𝜙(T, g::gas)
   S = 0
   for (key,val) in g.Y
      S = S + val * 𝜙(T, spd[key])
   end
   return S
end


"""
Returns Δs from the reference point defined at
Tref = 298.15 K
Pref = 101325 Pa

using the entropy complement function and
the entropy change due to pressure

`sp` can be either of type `species` or `gas`
"""
function s(T, P, sp)
   Tref = 298.15
   Pref = 101325

   Δs = 𝜙(T, sp) - 𝜙(Tref, sp) - ℜ*log(P/Pref)

end



# Specific functions for gas Compression
PR = 10
p2, T2 = 101325, 298.15
ηp = 0.90

""" 
Adiabatic compression given the 
compression pressure ratio (`PR`), the initial pressure (`p`)
and initial temperature (`T`).

Returns `Tfinal` and `pfinal`
"""
function compress(PR, p, T)
   Tfinal = T * PR^(ℜ/cp(T,Air))

   for i in 1:10
      Res  = (𝜙(Tfinal, Air) - 𝜙(T, Air))/ℜ - log(PR)
      Res′ = cp(Tfinal,Air)/ℜ/Tfinal
      dT  = Res/Res′
      Tfinal = Tfinal - dT
      # println(Tfinal)
      if abs(dT) < ϵ
         break
      end
   end

   return Tfinal, p*PR

end
"""
Adiabatic with NL solve
i.e. find x such that F(x)=0
"""
T = 298.15
p = 101325.
PR = 2.0
function f(x)
   s(T,p,Air) - s(x[1],p*PR,Air)
end

"""
Compression with polytropic efficiency
"""
function compress(PR, p, T, ηp )
   Tfinal = T * PR^(ℜ/cp(T,Air)/ηp)

   for i in 1:25
      Res  = (𝜙(Tfinal, Air) - 𝜙(T, Air))/ℜ - log(PR)/ηp
      Res′ = cp(Tfinal,Air)/ℜ/Tfinal
      dT  = Res/Res′
      
      # ω = 1.0
      # if i>10
      #    ω = 0.5
      # end
      # Tfinal = Tfinal - ω*dT
      Tfinal = Tfinal - dT
      # println("$i: $Tfinal $dT")
      if abs(dT) < ϵ
         break
      end
   end

   return Tfinal, p*PR

end

using PyCall
pygui()
using PyPlot
pygui(true)
plt.style.use(["~/prash.mplstyle", "seaborn-colorblind"])

function plot_polycomp()
   p2 = 101.325e3
   T2 = 288.15
   fig, ax = plt.subplots(figsize = (8,5), dpi = 200)
   PR = LinRange(20,50,10)
   for eff in [0.85, 0.90, 0.91, 0.95]
      T = [T for (T,p) in compress.(PR, p2, T2, eff)]
      l, = ax.plot(PR, T); label = "\$ \\eta_p = $(eff*100)\$%"
      ax.text(PR[end], T[end], label, color = l.get_color())
   end
   
   model_point = ax.scatter(32.64, 827.25,label = "MIT CFM56-5B NPSS model",
   color = "k", marker="*", s=80, zorder=10.1)

   GA_model = ax.scatter(28.6, 798.2, label = "GT CFM56-7B model", color = "k",
                      marker = ".", s = 50, zorder = 10.1)
   
   Engs = ["CFM56-7B27","CFM56-5B3", "LEAP-1B28", "PW1133G"]
   PR = [29.0,32.6, 42.0, 38.67]
   annotate_props = Dict("xycoords"=>("data", "axes fraction"), 
                      "arrowprops"=>Dict("arrowstyle"=>"->","connectionstyle"=>"arc3"), 
                      "ha"=>"center", "va" => "center","fontsize" => 8)
                     #  "bbox"=Dict("fc"=>"w", "ec"=>"none")

   # for (pr,eng) in zip(PR, Engs)
   #    ax.annotate(eng, xy=(pr,0), xytext=(pr,0.08),  zorder = 10)
   #    ax.axvline(pr, color = "k", lw = 1.0, alpha = 0.3, ls = "--", zorder = 9)
   # end
   ax.legend()
   ax.set_ylim(691.2326281862709, 1091.441513548923)

   ax.set_xlabel("\$\\pi_{\\mathrm{oo}}\$")
   ax.set_ylabel("T\$_{t3}\$[K]")
   ax.set_title("\$T_{t3}\$ vs \$\\pi_{oo}\$ at SLS conditions (\$T_{amb}\$ = 288.15 K)")
   plt.tight_layout()
end