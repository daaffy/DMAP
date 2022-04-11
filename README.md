# DMAP

<!-- [![Build Status](https://github.com/daaffy/DMAP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/daaffy/DMAP.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/daaffy/DMAP.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/daaffy/DMAP.jl) -->

[Diffusion map](https://www.sciencedirect.com/science/article/pii/S1063520306000546) approach to the approximation of the [dynamic Laplacian](https://arxiv.org/abs/1411.7186). 

## Instructions
To use the project, clone the repository to your own machine:
```bash
git clone https://github.com/daaffy/DMAP
cd ./DMAP
```
The project has multiple dependendencies (e.g., LinearAlgebra.jl, DifferentialEquations.jl) indicated in the Project.toml file. To set up the environment from the Julia REPL simply call
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
From there you can go ahead and run any of the given examples.
