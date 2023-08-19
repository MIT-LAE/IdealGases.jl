
"""
Function to create the required temperature array
"""
function Tarray(T)
   return [T^-2, T^-1, 1.0, T, T^2, T^3, T^4, log(T)]
end


"""
    Tarray!(T, TT)

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
    thermo_table(gas::Gas, 
    Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0,
    print::Bool=true)

Quickly generate a table of cp, h and s for a gas
"""
function thermo_table(gas::Gas; 
    Tstart::Float64=Tstd, Tend::Float64=2000.0, Tinterval::Float64=100.0,
    print::Bool=false)
    
   Trange = range(Tstart, Tend, step=Tinterval)

   cp_array = zero(Trange)
   h_array  = zero(Trange)
   ϕ_array  = zero(Trange)
   if print
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
   else
      for (i,T) in enumerate(Trange)
         gas.T = T
         cp_array[i]= gas.cp
         h_array[i] =  gas.h
         ϕ_array[i] =  gas.𝜙
      end
   return Trange, cp_array, h_array, ϕ_array
   end
end