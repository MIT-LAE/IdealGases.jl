using Pkg
Pkg.activate("..")
using IdealGases

push!(LOAD_PATH, "../src")

using Documenter

makedocs(
    remotes = nothing,
    sitename = "IdealGases.jl documentation",
    pages = ["Home" =>"index.md",
    "Thermodynamic Data" => "readthermo.md",
    "Gas type" => "gas.md",
    "Ideal gas thermo" => "idealgasthermo.md",
    "Combustion"=>"combustion.md"]

)