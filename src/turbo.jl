

"""
    PressureRatio(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

Generic pressure ratio conversion. See [`compress`](@ref) and [`expand`](@ref)
"""
function PressureRatio(gas::AbstractGas, PR::Float64, ηp::Float64=1.0,)

    T0 = gas.T
    ϕ0 = gas.ϕ
    P0 = gas.P
    R = gas.R
 
    Tfinal = T0 * PR^(Runiv/gas.cp/ηp)
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