
"""
    fuelbreakdown(fuel::String)

Returns the number of C, H, O, and N atoms that the fuel is composed of.
"""
function fuelbreakdown(fuel::String)
    C,H,O,N = 0, 0, 0, 0
    chunks = [fuel[idx] for idx in findall(r"[A-Z][a-z]?\d*", fuel)]
    for chunk in chunks
        element, number = match(r"([A-Z][a-z]?)(\d*)", chunk).captures

        if isempty(number)
            number = 1
        else
            number = parse(Int, number)
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
CᵢHⱼOₖNₗ + n(O2) * O2 ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2
⟹ CᵢHⱼOₖNₗ           ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2 - n(O2) * O2 
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
    ΔY = [nCO2,   nN2,  nH2O,  nO2]
    # ΔY = [M_CO2, M_N2, M_H2O, M_O2] .* 
    #      [nCO2,   nN2,  nH2O,  nO2] ./ Mfuel

    # return as dict to make it easier to set Gas mass fractions
    names = ["CO2", "N2", "H2O", "O2"]
    Ydict = Dict(zip(names, ΔY))

    return Ydict

end

