# Introduction

!!! warning "JMD.jl is still in the early development phase"
    This package is still very early in its development, and there are more mature molecular dynamics packages in Julia. For instance, [Molly.jl](https://juliamolsim.github.io/Molly.jl/stable/) and [NQCDynamics.jl](https://nqcd.github.io/NQCDynamics.jl/stable/) both offer molecular dyanmics in Julia.

Julia Molecular Dynamics (`JMD.jl`) aims to offer users a simple, intuitive and easy-to-use molecular dynamics enviornment. It draws inspiration from Python's [ASE](https://wiki.fysik.dtu.dk/ase/index.html), but is intended to be faster and offer users more flexibility. The flexibility comes from the relative ease with which users can add their own methods to dynamics or other components of `JMD.jl`.

### Installation

`JMD.jl` is not yet on the general registry, so for now installation can be done via GitHub.

```julia-repl
pkg> add https://github.com/Cavenfish/JMD
```

If you are more adventerous, you can consider installing the `dev` branch of `JMD.jl`. This will get updates more frequently, which gives users more features but also comes with increased chances of bugs. 

```julia-repl
pkg> add https://github.com/Cavenfish/JMD#dev
```

### Features

Currently, `JMD.jl` is able to perform the following simulations/calculations on molecular systems.
  
  - Geometry optimizations
  - Harmonic frequency calculations
  - Classical molecular dynamics in the NVE and NVT ensemble