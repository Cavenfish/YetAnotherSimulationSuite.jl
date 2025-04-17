# Custom Potentials

Using a custom potential within JMD is fairly straightforward. A custom potential with full functionality requires you make 5 functions, however, depending on your usecase you can make fewer. 

### Initializer Function

This function needs to initialize the variables used within your custom potential. The function should take a single parameter, either a vector of Atoms or a JMD cell object, and return a struct containing the potential's variables.

```julia
using JMD

# Your Potential's Parameters
struct _EXP_PotVars{F <: Float64} <: JMD.PotVars
  a::F
  b::F
  c::F
end

# Initializer Function for Vector of Atoms
EXP(bdys::Vector{JMD.MyAtoms}) = _EXP_PotVars(
  1.2,
  2.2,
  3.2
)

# Initializer Function for JMD Cell 
EXP(cell::JMD.MyCell) = _EXP_PotVars(
  1.2,
  2.2,
  3.2
)
```

### Dynamics Function

```julia
function EXP(dv, v, u, p, t)

  # Loop over all molecules
  for mol in p.mols
    # Your custom intramolecular potential
  end

  # Loop over all pairs of molecules 
  for par in p.pars
    # Your custom intermolecular potential
  end

  # Update the accelerations 
  dv .= F ./ p.m

  # If running an NVT, apply the thermostat
  if p.NVT
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  # Store energy and forces at this timestep 
  # This is optional
  push!(p.energy, E)
  push!(p.forces, F)
end
```

### Optimizing Function

```julia
function EXP(F, G, x0, p)

  # Loop over all molecules
  for mol in p.mols
    # Your custom intramolecular potential
  end

  # Loop over all pairs of molecules 
  for par in p.pars
    # Your custom intermolecular potential
  end

  # If gradient-based optimization, update gradient
  if G != nothing
    tmp = [j for i in forces for j in i]
    for i in 1:length(G)
      G[i] = -tmp[i]
    end
  end

  # If energy-based optimization, return the energy
  if F != nothing
    return E
  end
end
```