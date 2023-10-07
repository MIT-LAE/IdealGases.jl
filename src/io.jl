"""
    Base.show(io::IO, sp::species)

Shows a simplified representation of the `species` instance.
"""
function Base.show(io::IO, sp::species)
    println(io, "Species \"$(sp.name)\"\nMW = $(sp.MW) g/mol")
end
"""
    Base.show(io::IO, sp::composite_species)

Shows a simplified representation of the `composite_species` instance.
"""
function Base.show(io::IO, sp::composite_species)
    divider = "-"^(8*2+2*1)
    print(io, "Composite Species: \"$(sp.name)\"\nMW = $(sp.MW) g/mol\nwith composition:\n")
    @printf(io, "%8s  %8s\n", "Species", "Xᵢ")
    S = 0.0
    for (name, Xi) in sp.composition
        if Xi != 0
            @printf(io, "%8s  %8.5f\n", name, Xi)
            S += Xi
        end
    end
    println(divider)
    @printf(io, "%8s  %8.5f\n", "Σ", S)
end
"""
    Base.show(io::IO, gas::Gas)

Shows a simplified representation of the `Gas` instance.
"""
function Base.show(io::IO, gas::Gas)
    print(io, "Gas(T = $(gas.T) K; P = $(gas.P/1000.0) kPa; MW = $(gas.MW) g/mol)")
end

"""
    Base.print(io::IO, gas::Gas)

Pretty print for `Gas` instance
"""
function Base.print(io::IO, gas::Gas)
   @printf(io, "Ideal Gas at\n%3s = %8.3f K\n%3s = %8.3f kPa\n%3s = %8.3f J/K/kg\n%3s = %8.3f kJ/kg\n%3s = %8.3f kJ/K/kg",
     "T", gas.T, "P", gas.P/1000.0, "cp", gas.cp, "h", gas.h/1000.0, "s", gas.s/1000.0)
   println(io, "\n\nwith composition:")
   composition(gas,io)
end

"""
    composition(gas::Gas, io::IO=stdout)

Prints out the composition (Y_i) of the gas
"""
function composition(gas::Gas, io::IO=stdout)
   divider = "-"^(8*3+2*3+9)
   @printf(io, "%s\n", divider)
   @printf(io, "%8s  %8s  %8s  %9s\n", "Species", "Xᵢ", "Yᵢ", "MW[g/mol]")
   @printf(io, "%s\n", divider)
   for (name, Xi, Yi, mw) in zip(spdict.name, gas.X, gas.Y, spdict.MW)
      if Yi != 0
         @printf(io, "%8s  %8.3f  %8.3f  %9.3f\n", name, Xi, Yi, mw)
      end
   end
   @printf(io, "%s\n", divider)
   @printf(io, "%8s  %8.3f  %8.3f  %9.3f\n", "Σ", 
   sum(gas.X), sum(gas.Y), gas.MW)
end

"""
    print_thermo_table(gas::Gas; 
    Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0,)

# Examples
```julia-repl
julia> print_thermo_table(gas)
-----------------------------
 Species        Yᵢ  MW[g/mol]
-----------------------------
      N2     1.000     28.013
-----------------------------
     ΣYᵢ     1.000     28.013
 
   i     T[K]   cp[J/K/kg]     h[kJ/kg]   𝜙[kJ/K/kg]   s[kJ/K/kg]
----------------------------------------------------------------
   1   298.15    1039.6566       0.0000       6.8399       6.8399
   2   398.15    1043.9513     104.1279       7.1411       7.1411
   .    ...        ...          ...            ...          ...
  17  1898.15    1277.6602    1873.4391       8.9314       8.9314
  18  1998.15    1283.9203    2001.5241       8.9972       8.9972
```
"""
function print_thermo_table(gas::Gas; 
    Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0,
    massbasis::Bool=true)
    Trange = range(Tstart, Tend, step=Tinterval)
    print_thermo_table(gas, Trange, massbasis=massbasis)

end

"""
    print_thermo_table(gas::Gas, Trange::AbstractVector; massbasis::Bool=true)

method to print thermo table for a given range of temperatures `Trange` and 
allows one to specify whether output is desired on a mass or molar basis. 
"""
function print_thermo_table(gas::Gas, Trange::AbstractVector; massbasis::Bool=true)

    Trange, cp_array, h_array, 𝜙_array, s_array = thermo_table(gas, Trange)
    k = massbasis ? 1 : gas.MW/1000.0
    composition(gas)
    println(" ")
    divider = "-"^(4+8+12*4+4)
    if massbasis
        @printf("%4s %8s %12s %12s %12s %12s\n",
        "i",  "T[K]", "cp[J/K/kg]", "h[kJ/kg]", "𝜙[kJ/K/kg]", "s[kJ/K/kg]")
    else  
        @printf("%4s %8s %12s %12s %12s %12s\n",
        "i",  "T[K]", "cp[J/K/mol]", "h[kJ/mol]", "𝜙[kJ/K/mol]", "s[kJ/K/mol]")
    end

    println(divider)
    for (i,T) in enumerate(Trange)
        @printf("%4d %8.2f %12.4f %12.4f %12.4f %12.4f\n",
        i,  T, k*cp_array[i], k*h_array[i]/1000.0, 
        k*𝜙_array[i]/1000.0, k*s_array[i]/1000.0)
    end
end