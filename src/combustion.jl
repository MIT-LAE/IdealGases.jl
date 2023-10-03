
"""
    fuelbreakdown(fuel::String)

Returns the number of C, H, O, and N atoms that the fuel is composed of.
"""
function fuelbreakdown(fuel::String)
    C,H,O,N = 0.0, 0.0, 0.0, 0.0
    if !isempty(findall(r"[^cChHoOnN.^[0-9]",fuel))
        error("The input fuel string $fuel contains
        elements other than C,H,O, and N.")
        return nothing
    end
    chunks = [fuel[idx] for idx in findall(r"[a-zA-Z][a-z]?\d*\.?\d*", fuel)]
    for chunk in chunks
        element, number = match(r"([a-zA-Z][a-z]?)(\d*\.?\d*)", chunk).captures
        element = uppercase(element)
        if isempty(number)
            number = 1
        else
            number = parse(Float64, number)
        end
        if element == "C"
            C = C + number
        elseif element == "H"
            H = H + number
        elseif element == "O"
            O = O + number
        elseif element == "N"
            N = N + number
        else
            error("Fuel can only contain C, H, O or N atoms!")
        end
    end
    return([C, H, O, N])

end 

"""
    reaction_change_fraction(fuel::String)

Returns the mass fraction change due to complete combustion

Assume fuel of type  CᵢHⱼOₖNₗ , then
```
  CᵢHⱼOₖNₗ + n(O2) * O2 ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2
⟹CᵢHⱼOₖNₗ              ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2 - n(O2)*O2 
```
# Examples
```julia-repl
julia> reaction_change_fraction("CH4")
Dict{String, Float64} with 4 entries:
  "O2"  => -3.98926
  "H2O" => 2.24595
  "CO2" => 2.74331
  "N2"  => 0.0
```
"""
function reaction_change_fraction(fuel::String)
    CHON = fuelbreakdown(fuel) # Returns number of C, H, O, and N atoms in fuel
    element_mass_array = [12.01070, 1.00794, 15.99940, 14.00670]

    Mfuel = dot(CHON, element_mass_array)

    M_CO2 = 44.00950
    M_N2  = 28.01340
    M_H2O = 18.01528
    M_O2  = 31.99880

    # Calculate the number of moles of products + O₂
    nCO2 = CHON[1]*1.0
    nN2  = CHON[4]/2
    nH2O = CHON[2]/2
    nO2  = - (CHON[1] + CHON[2]/4 - CHON[3]/2) # Oxygen is used up/ lost

    ΔY = [M_CO2, M_N2, M_H2O, M_O2] .* 
         [nCO2,   nN2,  nH2O,  nO2] ./ Mfuel

    # return as dict to make it easier to set Gas mass fractions
    names = ["CO2", "N2", "H2O", "O2"]
    Ydict = Dict(zip(names, ΔY))

    return Ydict

end

"""
"""
function reaction_change_molar_fraction(fuel::AbstractString)
    CHON = fuelbreakdown(fuel) # Returns number of C, H, O, and N atoms in fuel
    # Calculate the number of moles of products + O₂
    nCO2 = CHON[1]*1.0
    nN2  = CHON[4]/2
    nH2O = CHON[2]/2
    nO2  = - (CHON[1] + CHON[2]/4 - CHON[3]/2) # Oxygen is used up/ lost

    X = [nCO2, nN2, nH2O, nO2]
    # return as dict to make it easier to set Gas mass fractions
    # names = ["CO2", "N2", "H2O", "O2"]
    # Xdict = Dict(zip(names, X))
    return X
end  # function reaction_change_molar_fraction


"""
"""
function vitiated_mixture(fuel::species, oxidizer::species, 
    FAR::Float64, ηburn::Float64=1.0)

    molFAR = FAR * oxidizer.MW/fuel.MW

    nCO2, nN2, nH2O, nO2 = ηburn .* molFAR .* reaction_change_molar_fraction(fuel.name)

    names = [oxidizer.name,   fuel.name, "CO2", "H2O", "N2", "O2"]
    ΔX    = [          1.0, 1.0 - ηburn,  nCO2,  nH2O,  nN2,  nO2]

    Xdict = Dict(zip(names, ΔX))
    return Xdict
end  # function vitiated_gas

"""
"""
function vitiated_species(fuel::species, oxidizer::species, 
    FAR::Float64, ηburn::Float64=1.0, name::AbstractString="vitiated species")
    
    names = spdict.name
    Xdict = vitiated_mixture(fuel, oxidizer, FAR, ηburn)
    X = zeros(MVector{length(names), Float64})
    for (key,value) in Xdict
       index = findfirst(x->x==key, names)
       X[index] = value
    end

    return generate_composite_species(X,name)
end  # function vitiated_species

function vitiated_species(fuel::AbstractString, oxidizer::AbstractString, 
    FAR::Float64, ηburn::Float64=1.0, name::AbstractString="vitiated species")

    o = spdict[findfirst(x->x==oxidizer, spdict.name)]
    f = spdict[findfirst(x->x==fuel, spdict.name)]
    return vitiated_species(f, o, 
    FAR, ηburn,  name)

end