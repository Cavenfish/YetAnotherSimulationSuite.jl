# Custom Potentials

Using a custom potential within `JMD.jl` is fairly straightforward. You need to make the following components.

  - `PotVars` struct
  - Initializer function (minimum 1, maximum 2)
  - Evaluation functions (minimum 1, maximum 3)

### Potential Variables Struct

This struct should hold all parameters (outside of particle positions, velocities, etc.) needed to evaluate your potential. The name of this struct is not important but the type must be the `PotVars` type from `JMD`. 

```julia
using JMD 

struct MyExamplePotVars{F<:Float64} <: JMD.PotVars
  gamma::F
end
```

### Initializer Function

This function needs to initialize the variables used within your custom potential. The function should take a single parameter, either a vector of `Particle` or a `JMD cell` object, and return your custom `PotVars` struct initialized with your potential's parameters.

```julia
MyExamplePotential(x::Union{Vector{JMD.MyAtoms}, JMD.MyCell}) = MyExamplePotVars(16.0)
```

### Evaluation functions

`JMD.jl` calculators can have 3 different evaluation functions: 
  
  - Energy only evaluation
  - Energy and non-inplace forces evaluation
  - Energy and inplace forces evaluation

