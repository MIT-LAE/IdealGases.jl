# Turbomachinery related functions

"""
    PressureRatio(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

Generic pressure ratio conversion. See [`compress`](@ref) and [`expand`](@ref)
"""
function PressureRatio(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

    T0 = gas.T
    ϕ0 = gas.ϕ
    P0 = gas.P
    R = gas.R
 
    Tfinal = T0 * PR^(R/gas.cp/ηp)
    Pfinal = P0*PR
    dT = Tfinal
    gas.P = Pfinal
    gas.T = Tfinal
    logPR_ηp = log(PR)/ηp
 
    for i in 1:20# abs(dT)>ϵ
       res  = (gas.ϕ - ϕ0)/R - logPR_ηp 
       res_dT = gas.ϕ_T/R
       dT  = - res/res_dT
 
       if abs(dT) ≤ ϵ
          break
       end
 
       Tfinal = Tfinal + dT
       gas.T = Tfinal
       # println("$i: $Tfinal $dT")
    end
 
    if abs(dT) > ϵ
       error("Error: PressureRatio did not converge:\ngas=", print(gas),
        "\n\nabs(dT) = ", abs(dT), " > ϵ (", ϵ, ")")
    end
 
    return gas
 
 end
 
 """
     compress(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)
 
 Compression with an optional polytropic efficiency.PR should be ≥ 1.0.
 See also [`expand`](@ref).
 """
 function compress(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)
    if PR < 1.0
       error("The specified pressure ratio (PR) to compress by needs to be ≥ 1.0.
      Provided PR = $PR. Did you mean to use `expand`?")
    end
    return PressureRatio(gas, PR, ηp)
 end
 
 """
     expand(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)
 
 Expansion at a given polytropic efficiency. PR should be ≤ 1.0.
 See also [`compress`](@ref).
 """
 function expand(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)
    if PR>1.0
       error("The specified pressure ratio (PR) to compress by needs to be ≤ 1.0.
      Provided PR = $PR. Did you mean to use `compress`?")
    end
    return PressureRatio(gas, PR, 1/ηp)
 end

"""
      gas_Mach!(gas::AbstractGas, M0::Float64, M::Float64, ηp::Float64 = 1.0)

Calculates the gas state for a change in Mach number with an optional polytropic efficiency.
 """
 function gas_Mach!(gas::AbstractGas, M0::Float64, M::Float64, ηp::Float64 = 1.0)
   itmax = 10
   ttol = 0.000001

   uosq = M0^2 * gas.γ * gas.R * gas.T
   h0 = gas.h
   R0 = gas.R 
   ϕ0 = gas.ϕ

   #---- initial guess for temperature, using constant-gamma relation
   t = gas.T * (1.0 + 0.5 * (gas.γ - 1) * M0^2) /
       (1.0 + 0.5 * (gas.γ - 1) * M^2)

   dt = 0.0
   #---- Newton iteration for actual temperature
   for iter = 1:itmax
      set_TP!(gas, t, gas.P) #keep pressure constant for now; will be changed later
      
      #------ usq( t, cp(t,al), r[al] , m )
      usq = M^2 * gas.γ * gas.R * gas.T
      usq_t = M^2 * gas.γ * gas.R
      usq_cp = M^2 * gas.R / (gas.cp - gas.R) * gas.T - usq / (gas.cp - gas.R)

      #------ res( h[t], usq( t, cp[t] ) )
      res = gas.h + 0.5 * usq - h0 - 0.5 * uosq
      res_t = gas.cp + 0.5 * (usq_t + usq_cp * gas.cp_T)
      #
      dt = -res / res_t

      if (abs(dt) < ttol)
         P = gas.P * exp(ηp * (gas.ϕ - ϕ0) / R0) #calculate pressure with polytropic efficiency
         set_TP!(gas, gas.T, P)
         return
      end

      t = t + dt
   end
end