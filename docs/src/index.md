# Introduction

!!! warning "YetAnotherSimulationSuite.jl is still in the pre-release phase"
    This package is still very early in its development, and there are more mature molecular dynamics packages in Julia. For instance, [Molly.jl](https://juliamolsim.github.io/Molly.jl/stable/) and [NQCDynamics.jl](https://nqcd.github.io/NQCDynamics.jl/stable/) both offer molecular dyanmics in Julia.

`YetAnotherSimulationSuite.jl` (YASS) aims to offer users a simple, intuitive and easy-to-use molecular dynamics enviornment. It draws inspiration from Python's [ASE](https://wiki.fysik.dtu.dk/ase/index.html), but is intended to be faster and offer users more flexibility. The flexibility comes from the relative ease with which users can add their own methods to dynamics or other components of YASS.

### Installation

YASS can be installed using the Julia package manger via:

```julia-repl
pkg> add YetAnotherSimulationSuite
```

If you are more adventerous, you can consider installing YASS from GitHub. This will get updates more frequently, which gives users more features but also comes with increased chances of bugs. 

```julia-repl
pkg> add https://github.com/Cavenfish/YetAnotherSimulationSuite.jl
```

### Features

Currently, YASS is able to perform the following simulations/calculations on molecular systems.
  
  - Geometry optimizations
  - Harmonic frequency calculations
  - Classical molecular dynamics in the NVE and NVT ensemble