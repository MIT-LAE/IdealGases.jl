struct species
    a_dict
    MW::Float64
    Hf::Float64
end

"""
Reads a NASA 9 polynomial thermo definintion file which can be obtained from
[NASA thermobuild](https://cearun.grc.nasa.gov/ThermoBuild/index_ds.html) 
and returns a dictionary of species
"""
function readThermo(filename)

    f = open(filename)
    input = readlines(f)

    # Get temperatures
    Temps = parse.(Float64, filter!(!isempty, split(strip(input[2]), " "))[1:end-1])[1:end-1]
    Tintervals = length(Temps)
    
    # Count number of species names and lines
    sp_names, sp_lines = zip([[strip(i[1:13]),n] for (n,i) in enumerate(input) if startswith(i, r"[a-zA-Z]\w[^END][^thermo]") ]...)
    Species = Dict()
    Nspecies = length(sp_names)

    a = zeros((2, 9)) #Assume just 2 temp intervals for now
    for (i, sp) in enumerate(sp_names)
        istart = sp_lines[i]
        MW = parse(Float64, strip(input[istart+1][53:65]))
        Hf = parse(Float64, strip(input[istart+1][66:80]))
        
        # Create arraay to store sp_names

        
        coeffs = replace(input[istart+3], "D"=>"e")*replace(input[istart+4], "D"=>"e")
        a[1,:]  = parse.(Float64, [coeffs[(i-1)*16+1:i*16] for i in range(1,10, step = 1) if i != 8])

        coeffs = replace(input[istart+6], "D"=>"e")*replace(input[istart+7], "D"=>"e")
        a[2,:]  = parse.(Float64, [coeffs[(i-1)*16+1:i*16] for i in range(1,10, step = 1) if i != 8])

        Species[sp] = species(Dict(200.0=>a[1,:], 1000.0=>a[2,:]), MW, Hf)
    end
    return Species
end
