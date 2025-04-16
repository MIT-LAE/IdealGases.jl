# IdealGasThermo.jl

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://mit-lae.github.io/IdealGasThermo.jl/dev/) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mit-lae.github.io/IdealGasThermo.jl/stable/) ![Codecov](https://img.shields.io/codecov/c/github/MIT-LAE/IdealGasThermo.jl) [![CI](https://github.com/MIT-LAE/IdealGasThermo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MIT-LAE/IdealGasThermo.jl/actions/workflows/CI.yml)


A simple package to efficiently calculate ideal gas properties based on NASA polynomials. This assumes that the specifc heat
of the gas/mixture is only a function of its temperature, i.e., $c_p(T)$ , $h(T)$, and $s(T,p)$ (note the entropy is a function of both pressure and temperature).

This package also introduces a method for fast calculations of thermodynamic properties of mixtures with fixed composition by calculating an equivalent polynomial representation of the mixture.
