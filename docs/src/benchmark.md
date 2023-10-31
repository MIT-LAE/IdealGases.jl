## [Performance](@id performance)

In general calculating the thermodynamic properties of a mixture scales with the number of constituent species. 
The key feature of `IdealGases.jl` is the ability to represent
a mixture of gases as a single `composite_species` which can provide noticable performance improvements.

Below is an example. We first create a function that sets the gas temp 100 times and calculates $c_p$, $h$, $\phi$, and $dc_p/dT$. We do this so that 


```julia
TT = rand(200.:600., 100)

function benchmark_Gas(TT::AbstractVector, gas::AbstractGas)
    @views for i in eachindex(TT)
        gas.T = TT[i]; 
        gas.cp
        gas.ϕ
        gas.h
        gas.cp_T
    end
end
```

We now initialze a `Gas` and `Gas1D` instances and benchmark the above function.

```julia
gas = Gas()
gas1D = Gas1D()

using BenchmarkTools

julia> @benchmark benchmark_Gas($TT, $gas1D)
BenchmarkTools.Trial: 10000 samples with 9 evaluations.
 Range (min … max):  2.778 μs …  12.551 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.884 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.909 μs ± 202.135 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▂█▅▅ ▄▇▇▅▃▃▃▃▂▂▁                                            ▂
  ████████████████▆▆▇▇▆▆▆██▆▆▇▇█▇▇▇▇▇▇▆▇▆▆▅▆▆▅▅▆▄▆▄▆▅▆▅▆▆▆▅▄▆ █
  2.78 μs      Histogram: log(frequency) by time      3.74 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark benchmark_Gas($TT, $gas)
BenchmarkTools.Trial: 10000 samples with 4 evaluations.
 Range (min … max):  7.604 μs …  13.541 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     7.667 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   7.922 μs ± 561.794 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ██▂ ▂▄▅▄▅▃▁▂▂▂                                              ▂
  ███▆███████████▇▇▇█▇▇▇▇█▇▆██▇███▇▆▇▆▇▇▇▆▇▇▇▆▇▆▇▆▆▇▇▆▆▇▇▇▆▅▅ █
  7.6 μs       Histogram: log(frequency) by time      10.2 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

```

We see that the `Gas1D` calculations are ~3x faster. If we repeat this with a mixture that has more species the time taken for the `Gas` type calculation roughly scales with the number of species where as the `Gas1D` type calculations are independent of the number of constituent species in the gas.

For 10 components:

```julia-repl
julia> @benchmark benchmark_Gas($TT, $gas)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  17.208 μs … 149.167 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     17.875 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   18.134 μs ±   1.755 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▂▃▁  ▁█▇▅▃▂▃▄▄▃▂▁▁▁▂▃▁                                       ▂
  ███▆▅█████████████████▇▅▅▅▆▅▅▅▄▅▅▅▅▄▄▅▅▅▆▃▅▁▅▅▆▄▄▅▅▅▅▆▅▅▄▅▆▆ █
  17.2 μs       Histogram: log(frequency) by time      22.6 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark benchmark_Gas($TT, $gas1D)
BenchmarkTools.Trial: 10000 samples with 9 evaluations.
 Range (min … max):  2.792 μs …  13.907 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.903 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.958 μs ± 233.528 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▁▁   █▇▆▄▄▄▄▂▁▁▂▂▂▁                                         ▁
  ██▇▄███████████████▇▅▅▆▆▆▆▆▆▄▆▇▇▇▇▆▆▇▇▇▇▆▆▆▆▅▆▆▆▅▇▅▆▅▅▅▅▅▆▆ █
  2.79 μs      Histogram: log(frequency) by time      3.75 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.
```