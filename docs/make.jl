# using Pkg
# Pkg.activate("..")
using Documenter, IdealGases

push!(LOAD_PATH, "../src")


makedocs(
    remotes=nothing,
    sitename="IdealGases.jl",
    pages=["Home" => "index.md",
        "Thermodynamic Data" => "readthermo.md",
        "Gas and Species" => "gas.md",
        "Ideal gas thermodynamics" => "idealgasthermo.md",
        "Combustion" => "combustion.md",
        "Performance benchmarks" => "benchmark.md"], 

    format = Documenter.HTML(; mathengine=
        Documenter.KaTeX(
            Dict(:delimiters => [
                Dict(:left => raw"$",   :right => raw"$",   display => false),
                Dict(:left => raw"$$",  :right => raw"$$",  display => true),
                Dict(:left => raw"\[",  :right => raw"\]",  display => true),
                ],
                :macros => 
                Dict("\\RR" => "\\mathbb{R}",
                    "\\genfuel" => 
                    raw"{\mathrm{C}_{x_{\mathrm{C}}}\mathrm{H}_{x_{\mathrm{H}}}\mathrm{O}_{x_{\mathrm{O}}}\mathrm{N}_{x_{\mathrm{N}}}}",
                    raw"\fst" => raw"f_{\mathrm{stoich.}}",
                ),
            )
        )
    )
)
