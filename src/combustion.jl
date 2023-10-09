
"""
    fuelbreakdown(fuel::String)

Returns the number of C, H, O, and N atoms that the fuel is composed of.

# Examples
```julia-repl

julia> IdealGases.fuelbreakdown("CH4")' #transpose is simply to save space in the docs
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 1.0  4.0  0.0  0.0

julia> IdealGases.fuelbreakdown("C12H23.5")'
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 12.0  23.5  0.0  0.0

julia> IdealGases.fuelbreakdown("CH3COOH")'
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 2.0  4.0  2.0  0.0

julia> IdealGases.fuelbreakdown("CH3CH2OH")'
1×4 adjoint(::Vector{Float64}) with eltype Float64:
 2.0  6.0  1.0  0.0

```
"""
function fuelbreakdown(fuel::String)
    C,H,O,N = 0.0, 0.0, 0.0, 0.0
    if !isempty(findall(r"[^cChHoOnN.^[0-9]",fuel))
        try
            fuel = species_in_spdict(fuel).formula 
        catch e
            if isa(e, ArgumentError)
                error("""The input fuel string $fuel is not found in 
                the thermo database and contains
                elements other than C,H,O, and N.\n""")
            end
        end
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

fuelbreakdown(fuel::species) = fuelbreakdown(fuel.formula)

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
    reaction_change_molar_fraction(fuel::AbstractString)

Returns the mole fraction change due to complete combustion of one mole of
the specified fuel

Assume fuel of type  CᵢHⱼOₖNₗ , then
```
    CᵢHⱼOₖNₗ + n(O2) * O2 ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2
   ⟹CᵢHⱼOₖNₗ              ---> n(CO2)*CO2 + n(H2O)*H2O + n(N2)*N2 - n(O2)*O2 
```

# Examples
```julia-repl
julia> IdealGases.reaction_change_molar_fraction("CH4")
4-element Vector{Float64}:
  1.0
  0.0
  2.0
 -2.0
```
"""
function reaction_change_molar_fraction(fuel::AbstractString)
    CHON = fuelbreakdown(fuel) # Returns number of C, H, O, and N atoms in fuel
    # Calculate the number of moles of products + O₂
    nCO2 = CHON[1]*1.0
    nN2  = CHON[4]/2
    nH2O = CHON[2]/2
    nO2  = - (CHON[1] + CHON[2]/4 - CHON[3]/2) # Oxygen is used up/ lost

    X = [nCO2, nN2, nH2O, nO2]
    return X
end  # function reaction_change_molar_fraction

"""
    stoich_molar_fuel_oxy_ratio(fuel::AbstractString)

Calculates the molar fuel-**oxygen** ratio for stoichiometric combustion.

# Examples
```julia-repl
julia> using IdealGases

julia> IdealGases.stoich_molar_fuel_oxy_ratio("CH4")
0.5

julia> IdealGases.stoich_molar_fuel_oxy_ratio("C12H23")
0.056338028169014086
```
"""
function stoich_molar_fuel_oxy_ratio(fuel::AbstractString)
    X = reaction_change_molar_fraction(fuel)
    molFuelOxyRatio = abs(1/X[end])
    return molFuelOxyRatio
end  # function stoich_molar_FAR

"""
    stoich_molar_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies)

Calculates the **molar** fuel-oxidizer ratio for stoichiometeric combustion for 
and arbitrary fuel and oxidizer.

# Examples
```julia-repl
julia> CH4 = species_in_spdict("CH4");

julia> IdealGases.stoich_molar_FOR(CH4)
0.104738
```
"""
function stoich_molar_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies=DryAir)
    molFuelOxyRatio = stoich_molar_fuel_oxy_ratio(fuel.name)
    if typeof(oxidizer)== species
        Xin = Dict(oxidizer.name => 1.0)
        if oxidizer.name == "Air"
            Xin = Xair
        end
    else
        Xin = oxidizer.composition
    end

    return molFuelOxyRatio * Xin["O2"]

end  # function stoich_molar_FOR 

"""
    stoich_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies=DryAir)

Calculates the **mass** fuel-oxidizer ratio for stoichiometeric combustion for 
and arbitrary fuel and oxidizer.

# Examples
```julia-repl
julia> IdealGases.stoich_FOR(CH4)
0.05800961333050494
```
"""
function stoich_FOR(fuel::AbstractSpecies, oxidizer::AbstractSpecies=DryAir)
    molFOR = stoich_molar_FOR(fuel, oxidizer)
    return molFOR * fuel.MW/oxidizer.MW
end  # function stoich_FOR

stoich_FOR(fuel::AbstractString, oxi::AbstractString) = 
stoich_FOR(species_in_spdict(fuel), species_in_spdict(oxi))

"""
    vitiated_mixture(fuel::AbstractSpecies, oxidizer::AbstractSpecies, 
    FAR::Float64, ηburn::Float64=1.0)

Calculates the composition of a burnt gas mixture. Defaults to stoichiometric 
conditions if FAR is not specified. `vitiated_mixture` returns the number of 
moles of each species present in the burnt gas mixture after combustion at the 
specified FAR. *Note* the sum of result in general will **not** sum to 1. 

# Examples
```julia-repl
julia> CH4 = species_in_spdict("CH4");

julia> Air = IdealGases.DryAir;

julia> IdealGases.vitiated_mixture(CH4, Air, 0.04)
Dict{Any, Any} with 6 entries:
  "O2"  => 0.0650337
  "CH4" => 0.0
  "Ar"  => 0.009365
  "H2O" => 0.144442
  "CO2" => 0.0725401
  "N2"  => 0.78084

julia> IdealGases.vitiated_mixture(CH4, Air)
Dict{Any, Any} with 6 entries:
  "O2"  => 0.0
  "CH4" => 0.0
  "Ar"  => 0.009365
  "H2O" => 0.209476
  "CO2" => 0.105057
  "N2"  => 0.78084
```
See [here](@ref vitiated) for some explanation of the background.
"""
function vitiated_mixture(fuel::AbstractSpecies, oxidizer::AbstractSpecies, 
    FAR::Float64, ηburn::Float64=1.0)

    if typeof(oxidizer)== species
        Xin = Dict(oxidizer.name => 1.0)
        if oxidizer.name == "Air"
            Xin = Xair
        end
    else
        Xin = oxidizer.composition
    end

    molFAR = FAR * oxidizer.MW/fuel.MW

    nCO2, nN2, nH2O, nO2 = ηburn .* molFAR .* reaction_change_molar_fraction(fuel.name)

    nFuel = molFAR*(1.0 - ηburn)

    names = [fuel.name, "CO2", "H2O", "N2", "O2"]
    ΔX    = [    nFuel,  nCO2,  nH2O,  nN2,  nO2]

    Xdict = Dict(zip(names, ΔX))
    
    Xvitiated = mergewith(+, Xin, Xdict)

    return Xvitiated
end  # function vitiated_gas

vitiated_mixture(fuel, oxidizer) = 
vitiated_mixture(fuel, oxidizer, stoich_FOR(fuel, oxidizer))

"""
    vitiated_mixture(fuel::AbstractString, oxidizer::AbstractString, 
    FAR::Float64, ηburn::Float64=1.0)

Convinence function that finds fuel and oxidizer from thermo database
"""
function vitiated_mixture(fuel::AbstractString, oxidizer::AbstractString, 
    FAR::Float64, ηburn::Float64=1.0)

    fuel = species_in_spdict(fuel)
    oxidizer = species_in_spdict(oxidizer)

    return vitiated_mixture(fuel, oxidizer, FAR, ηburn)
    
end  # function vitiated_mixture

"""
    vitiated_species(fuel::AbstractSpecies, oxidizer::AbstractSpecies, 
    FAR::Float64, ηburn::Float64=1.0, name::AbstractString="vitiated species")

Returns a [`composite_species`](@ref) that represents the burnt gas mixture
at the specified FAR. If no FAR is provided stoichiometeric conditions are
assumed. 

# Examples
```julia-repl
julia> IdealGases.vitiated_species(CH4, Air, 0.05, name = "CH4-Air-0.05")
Composite Species: "CH4-Air-0.05"
MW = 27.89510190126262 g/mol
with composition:
 Species        Xᵢ
      O2   0.02653
      Ar   0.00859
     H2O   0.16560
     CO2   0.08309
      N2   0.71619
------------------
       Σ   1.00000
```
"""
function vitiated_species(fuel::AbstractSpecies, oxidizer::AbstractSpecies, 
    FAR::Float64; ηburn::Float64=1.0, name::AbstractString="vitiated species")
    
    Xdict = vitiated_mixture(fuel, oxidizer, FAR, ηburn)
    molFAR = FAR * oxidizer.MW/fuel.MW
    totalmoles = 1 + molFAR
    names = spdict.name

    X = zeros(MVector{length(names), Float64})
    for (key,value) in Xdict
       index = findfirst(x->x==key, names)
       X[index] = value/totalmoles
    end
    
    return generate_composite_species(X,name)
end  # function vitiated_species

function vitiated_species(fuel::AbstractString, oxidizer::AbstractString, 
    FAR::Float64; ηburn::Float64=1.0, name::AbstractString="vitiated species")

    fuel = species_in_spdict(fuel)
    oxidizer = species_in_spdict(oxidizer)
    return vitiated_species(fuel, oxidizer, FAR; ηburn=ηburn,  name=name)

end

vitiated_species(f::AbstractSpecies,o::AbstractSpecies) = 
vitiated_species(f, o, stoich_FOR(f, o))

vitiated_species(f::AbstractString,o::AbstractString) = 
vitiated_species(f, o, stoich_FOR(f, o))

"""
    fixed_fuel_vitiated_species(fuel, oxidizer, ηburn::Float64=1.0)

Returns a function `burntgas(FAR::Float64)` that is specific to the fuel and oxidizer
combination provided. This gives a highly performant function that can 
simply be called at any given FAR for that specific fuel+oxidizer combo.

This is ~4x faster than doing 
```burntgas(FAR) = vitiated_species("CH4", "Air", FAR)```

## Examples
```julia-repl
julia> burntgas = IdealGases.fixed_fuel_vitiated_species(CH4, Air)
(::IdealGases.var"#burntgas#52"{species, composite_species, Vector{Float64}, Vector{Float64}, Float64}) (generic function with 1 method)

julia> burntgas(0.05)
Composite Species: "burntgas(CH4 + Dry Air; 0.05)"
MW = 27.89510190126262 g/mol
with composition:
 Species        Xᵢ
      O2   0.02653
      Ar   0.00859
     H2O   0.16560
     CO2   0.08309
      N2   0.71619
------------------
       Σ   1.00000

julia> gas1 = Gas1D(burntgas(0.02)) #Returns a Gas1D intialized with the burnt gas properties
Gas1D(burntgas(CH4 + Dry Air; 0.02); MW = 28.51473501878705 g/mol)
at T = 298.15 K; P = 101.325 kPa
```
"""
function fixed_fuel_vitiated_species(fuel, oxidizer, ηburn::Float64=1.0)

    nCO2, nN2, nH2O, nO2 = ηburn .* reaction_change_molar_fraction(fuel.name)
    nFuel = (1.0 - ηburn)
    massratio =  oxidizer.MW/fuel.MW
    if typeof(oxidizer)== species
        Xin::Dict{String, Float64} = Dict(oxidizer.name => 1.0)
        if oxidizer.name == "Air"
            Xin = Xair
        end
    else
        Xin = oxidizer.composition
    end
    names::Vector{String} = [fuel.name, "CO2", "H2O", "N2", "O2"]
    ΔX::Vector{Float64}   = [nFuel,  nCO2,  nH2O,  nN2,  nO2]
    Xdict::Dict{String, Float64} = Dict(zip(names, ΔX))

    ΔX_array = zeros(Float64, Nspecies)
    Xin_array = zeros(Float64, Nspecies)
    allnames = view(spdict.name, :)

    for (key,value) in Xdict
        index::Int64 = findfirst(x->x==key, allnames)
        ΔX_array[index] = value
    end

    for (key,value) in Xin
        index::Int64 = findfirst(x->x==key, allnames)
        Xin_array[index] = value
    end

    """
       burntgas(FAR::Float64)
    Inner function that will be returned
    """
    function burntgas(FAR::Float64) 
        molFAR = FAR.*massratio
        X = ΔX_array.* molFAR
        @. X = (X + Xin_array)/(1+molFAR)
        generate_composite_species(X,
        "burntgas($(fuel.name) + $(oxidizer.name); $FAR)")

    end  # function burntgas

    return burntgas
end  # function fixed_fuel_vitiated_species
