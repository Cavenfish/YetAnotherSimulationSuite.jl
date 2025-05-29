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

  1) Energy only evaluation
  2) Inplace forces only evaluation
  3) Energy and inplace forces evaluation

You don't need to make all 3 functions, but you must have both energy and force evaluations to make a full calculator. This means you can make just function number 3 or both functions 1 and 2. If you make all functions the calculator is given more flexibility which may increase performance. For instance, if you want to get the potential energy of your system and only function number 3 is available then the calculator will need to evaluate both forces and energy to return just the energy. Whereas, an energy only evaluation may be faster and more memory efficient.

The evaluation functions have required arguments and not adhereing to them will cause the calculator to not work. In the examples below they are shown, their names however are allowed to changed.

```julia
"""
Energy only evaluation

Required arguments:
  u: positions of particles in system
  vars: struct of variables passed to the calculator

Required return:
  energy: the energy of the system

Order must be preserved but u and vars can be named anything.
"""
function MyEnergy(u, vars)
  # Initialize energy
  E = 0.0

  # Access your PotVars from vars
  potVars = vars.potVars

  # Vars also has the masses and species of the system
  m = vars.m # masses
  s = vars.s # species (ie. "H" for hydrogen)

  # Iterate over all molecules in your system
  for mol in vars.mols
    E += potVars.gamma
  end

  # Iterate over all pairs in your system
  for par in vars.pars
    E -= 1.0
  end

  # Return the energy
  E
end

"""
Inplace forces only evaluation

Required arguments:
  F: forces on particles in system
  u: positions of particles in system
  vars: struct of variables passed to the calculator

No return value

Order must be preserved but u and vars can be named anything.
"""
function MyForces!(F, u, vars)
  # This exmaple just illustrates the declration of
  # this function.
  F .= 0.0
end

"""
Energy and inplace forces only evaluation

Required arguments:
  F: forces on particles in system
  u: positions of particles in system
  vars: struct of variables passed to the calculator

Required return:
  energy: the energy of the system

Order must be preserved but u and vars can be named anything.
"""
function MyEnergyAndForces!(F, u, vars)
  # You are allowed to do this
  E = MyEnergy(u, vars)
  MyForces!(F, u, vars)

  # Keep in mind though that if you need to iterate
  # over all mols and pars for both the forces and energy
  # doing this only once in here would be faster and more
  # memory efficient. Here we call the other two function 
  # just for the sake of showing you can do such a thing.

  # Return energy
  E
end
```

### Putting it all Together

Now that we have made all the necessary components we can put it all together as a calculator.

```julia
# Define it as a variable
calc = Calculator(
  MyExamplePotential; 
  E = MyEnergy,
  F = MyForces!,
  EF = MyEnergyAndForces!
)

# Or define a grabber function for it
MyPotential() = Calculator(
  MyExamplePotential; 
  E = MyEnergy,
  F = MyForces!,
  EF = MyEnergyAndForces!
)

# You can then use either style
E = getPotEnergy(calc, bdys)
E = getPotEnergy(MyPotential(), bdys)
```
