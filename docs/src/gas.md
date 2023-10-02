# Representing ideal gases

`IdealGases.jl` exports the `Gas` type which stores relevant thermodynamic
information about the gas.

```@docs
Gas
Gas()
Base.setproperty!(gas::Gas, sym::Symbol, val::Float64)
```

## Setting the thermodynamic state of the gas

The following functions let you set the thermodynamic state of the gas. 
These functions change the state of the gas *in place* i.e., the gas object
is modified and no new copy is created.

```@docs
set_h!
set_hP!
set_TP!
set_Δh!
```