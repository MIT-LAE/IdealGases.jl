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
   @printf(io, "Ideal Gas at\n%3s = %8.3f K\n%3s = %8.3f kPa\n%3s = %8.3f J/K/mol\n%3s = %8.3f kJ/mol\n%3s = %8.3f kJ/K/mol",
     "T", gas.T, "P", gas.P/1000.0, "cp", gas.cp, "h", gas.h/1000.0, "s", gas.s/1000.0)
   println(io, "\n\nwith composition:")
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