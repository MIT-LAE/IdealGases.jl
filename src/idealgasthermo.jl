
"""
    Cp(Tarray::AbstractVector{T}, a::AbstractVector{T}) where T

Calculates cp of the given species in J/K/mol
(This is a completely non-allocating operation.)
```
Cp0/R = a₁T⁻² + a₂T⁻¹ + a₃ + a₄T + a₅T² + a₆T³ + a₇T⁴
```
"""
function Cp(Tarray::AbstractVector{T}, a::AbstractVector{T}) where T
    Cp_R = dot(view(a, 1:7), view(Tarray, 1:7))
    Cp = Cp_R*Runiv
    return Cp #J/K/mol
end

"""
    dCpdT(TT::AbstractVector{T}, a::AbstractVector{T}) where T

Returns the derivative dcp/dT [J/K²/mol]
```
dCp0/dT = R(-2a1*T^-3 -a2*T^-2 + a4 + 2a5*T + 3a6*T^2 + 4a7*T^3)
```
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
    h(TT::AbstractVector{type}, a::AbstractVector{type}) where type

Calculates h of the given **species** in J/mol

Calcualted by:
```
H0/RT = -a1*T^-2 + a2*T^-1*ln(T) + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + b1/T
      = -a1*T₁   + a2*T₂*T₈      + a3 + a4*T₄/2 + a5*T₅/3  + a6*T₆/4  + a7*T₇/5  + a₈*T₂
```
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
    𝜙(TT::AbstractVector{type},a::AbstractVector{type}) where type

Calculates the entropy complement function 𝜙=∫(cₚ/T)dT in J/K/mol

This is calculated at standard state. Tref = 298.15 K, Pref = 101325 Pa.
```
S0/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3.0 + a7*T^4/4 + b2 
     = -a1*T₁/2   - a2*T₂   + a3*T₈    + a4*T₄+ a5*T₅/2  + a6*T₆/3.0  + a7*T₇/4  + a₉   
```
"""
function 𝜙(TT::AbstractVector{type},a::AbstractVector{type}) where type
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

# For individual species:

"""
    Cp(T, sp::species)

Calculates cp for a **species** type in J/K/kg.
"""
function Cp(T, sp::species)
   TT = Tarray(T)
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   a = getfield(sp, s)
   return Cp(TT, a) * 1000.0/sp.MW

end

"""
    h(T, sp::species)

Calculates h for a species in J/kg
"""
function h(T, sp::species)
   TT = Tarray(T)
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   a = getfield(sp, s)
   h(TT, a) *1000.0/sp.MW
end

"""
    s(T, P, sp::species)

Calculates s for a species in J/K/kg
"""
function s(T, P, sp::species)
   TT = Tarray(T)
   if T<1000.0
      s = :alow
   else
      s = :ahigh
   end
   a = getfield(sp, s)
   sᵒ = 𝜙(TT, a) - Runiv*log(P/Pstd)
   return sᵒ * 1000.0/sp.MW
end