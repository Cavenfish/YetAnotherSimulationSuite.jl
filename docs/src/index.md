# Introduction

!!! warning "YASS.jl is still in the pre-release phase"
    This package is still very early in its development, and there are more mature molecular dynamics packages in Julia. For instance, [Molly.jl](https://juliamolsim.github.io/Molly.jl/stable/) and [NQCDynamics.jl](https://nqcd.github.io/NQCDynamics.jl/stable/) both offer molecular dyanmics in Julia.

Yet Another Simulation Suite (`YASS.jl`) aims to offer users a simple, intuitive and easy-to-use molecular dynamics enviornment. It draws inspiration from Python's [ASE](https://wiki.fysik.dtu.dk/ase/index.html), but is intended to be faster and offer users more flexibility. The flexibility comes from the relative ease with which users can add their own methods to dynamics or other components of `YASS.jl`.

### Installation

`YASS.jl` is not yet on the general registry, so for now installation can be done via GitHub.

```julia-repl
pkg> add https://github.com/Cavenfish/YASS
```

If you are more adventerous, you can consider installing the `dev` branch of `YASS.jl`. This will get updates more frequently, which gives users more features but also comes with increased chances of bugs. 

```julia-repl
pkg> add https://github.com/Cavenfish/YASS#dev
```

### Features

Currently, `YASS.jl` is able to perform the following simulations/calculations on molecular systems.
  
  - Geometry optimizations
  - Harmonic frequency calculations
  - Classical molecular dynamics in the NVE and NVT ensemble


### Dependencies

`YASS.jl` relies on several specialized external packages. These packages actively maintained and well trusted within the Julia ecosystem. If this changes, `YASS.jl` will remove these dependencies and if necessary implement the specialized code in-house. Here the dependencies are listed with links to their repos to give credit to their work but also to provide transperancy.

Packages within the Julia standard library are listed seperately since they are expected to be maintained as well as the Julia langague itself. For each dependency there is a short description of how it is used in `YASS.jl`, or why it is considered for removal.

**Julia Standard Library Packages**
  - [TOML](https://github.com/JuliaLang/julia/tree/master/stdlib/TOML):
  - [Libdl](https://github.com/JuliaLang/julia/tree/master/stdlib/Libdl):
  - [Statistics](https://github.com/JuliaLang/julia/tree/master/stdlib/Statistics):
  - [Serialization](https://github.com/JuliaLang/julia/tree/master/stdlib/Serialization):
  - [LinearAlgebra](https://github.com/JuliaLang/julia/tree/master/stdlib/LinearAlgebra):

**External Packages**
  - [FFTW](https://github.com/JuliaMath/FFTW.jl):
  - [Optim](https://github.com/JuliaNLSolvers/Optim.jl):
  - [PyCall](https://github.com/JuliaPy/PyCall.jl):
  - [Spglib](https://github.com/spglib/spglib):
  - [Chemfiles](https://github.com/chemfiles/chemfiles):
  - [Distances](https://github.com/JuliaStats/Distances.jl):
  - [Clustering](https://github.com/JuliaStats/Clustering.jl):
  - [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl):
  - [Distributions](https://github.com/JuliaStats/Distributions.jl):
  - [KernelDensity](https://github.com/JuliaStats/KernelDensity.jl):
  - [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl):

**Considered for Removal**
  - [JLD2](https://github.com/JuliaIO/JLD2.jl): This package is currently only used to load neural network data for potentials. This functionality can be covered by the `Serialization` package, which can reduce the total dependency count. Note, this is a well maintained package and users are encouraged to use it alongside `YASS.jl`.
  - [DataFrames](https://github.com/JuliaData/DataFrames.jl): This package does not add any functionality to `YASS.jl`, but rather enchances user exerpience. However, this can be achieved by users using the package alongside `YASS.jl`, rather than it being a dependency. Note, this is a well maintained package and users are encouraged to use it alongside `YASS.jl`.