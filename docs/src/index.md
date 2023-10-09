# IdealGases.jl Documentation

`IdealGases.jl` is a tool for fast thermodynamic properties calculations for
calorically/thermally perfect gases. This assumes that the specifc heat
of the gas/mixture is only a function of its temperature, i.e., $c_p(T)$ , 
$h(T)$, and $s(T,p)$.

One of the important features of `IdealGases.jl` is the ability to represent
a mixture of gases as a single [`composite_species`](@ref vitiated).


## Getting started

There are several workflows that are possible to use `IdealGases.jl`. We outline here the most common few.

### Simple install

The easiest way to run `IdealGases.jl` would be to add the package using the julia package manager using the github repository.

You can do this by starting a Julia session and then activating the package manager by typing `]` and then entering:
```julia-repl
pkg> add "git@github.mit.edu:LAE/Gas.jl.git"
```

You can then import `IdealGases` as you would with any Julia package:
```julia-repl
julia> using IdealGases
```
### Local development

If you are going to develop the source code of `IdealGases.jl` you might benefit from a local clone of the git repository which
can then fit into a workflow using [`Revise.jl`](https://timholy.github.io/Revise.jl/stable/) for example.

Step 1: Clone the git repo locally
```bash
git clone git@github.mit.edu:LAE/Gas.jl.git
```

Step 2: `cd` to the folder where IdealGases is cloned

Step 3: Use `Pkg` to install/ develop the package

```julia
pkg> dev .
```

You should now be able to import IdealGases from within any Julia script in your base environment.

!!! tip "Tip"
    If you are using `Revise.jl` be sure to first import `Revise` before importing `IdealGases`

    ```julia
    using Revise
    using IdealGases
    ```
