# Reading thermodynamic data

`IdealGasThermo.jl` uses the `NASA-9` polynomials to calculate the thermodynamic properties (``c_p, h, s``)

```@docs
readThermo
```
## Default species included

A number of species of interest are included in `IdealGasThermo.jl`.

```@example
using IdealGasThermo #hide
println(IdealGasThermo.spdict.name) #hide
```

## Adding new species

Any new gaseous species can be added for calculations by adding the NASA-9 polynomial data to the `thermo.inp` file. For example, if you wanted to add methane (this already exists in the database so you don't need to):
 - Step 1: Go to the [thermobuild](https://cearun.grc.nasa.gov/ThermoBuild/index_ds.html).
 - Step 2: Copy the data for the species of interest
   
```
   CH4               Gurvich,1991 pt1 p44 pt2 p36.                                 
 2 g 8/99 C   1.00H   4.00    0.00    0.00    0.00 0   16.0424600     -74600.000
    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10016.202
-1.766850998D+05 2.786181020D+03-1.202577850D+01 3.917619290D-02-3.619054430D-05
 2.026853043D-08-4.976705490D-12                -2.331314360D+04 8.904322750D+01
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10016.202
 3.730042760D+06-1.383501485D+04 2.049107091D+01-1.961974759D-03 4.727313040D-07
-3.728814690D-11 1.623737207D-15                 7.532066910D+04-1.219124889D+02
```

- Step 3: Append this to `thermo.inp`