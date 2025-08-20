# Thermostats

`YASS.jl` includes some pre-written thermostats for use in molecular dyunamics simulations. However, it is also fairly easy to create a custom thermostat to use in simulations. 

The included thermostats are:

  - Berendsen
  - Langevin
  - Canonical velocity rescaling


### Custom Thermostat

Here the `Berendsen` thermostat is shown to illustrate how to create a custom thermostat. Custom thermostats require a struct for their parameters and an action function.

```julia
struct Berendsen{F<:AbstractFloat} <:ThermoVars
  gamma::F
end

"""
Action function

Required arguments:
  a: accerleration of particles in system
  v: velocity of particles in system
  m: masses of particles in system
  Tsim: the current simulation temperature
  thermostat: the thermostat

return nothing

Order must be preserved but u and vars can be named anything.
"""
function Berendsen!(a, v, m, Tsim, thermostat)
  # Access our thermostat parameter gamma
  gamma = thermostat.vars.gamma

  # thermostat.T is the target temperature
  if Tsim == 0.0
    a .+= gamma .* v
  else
    a .+= gamma * (thermostat.T / Tsim - 1) .* v
  end
end
```

Now that we have made all the necessary components we can put it all together as a thermostat.

```julia
# Define a constructor for it passable temperature
function Berendsen(T::F, gamma::F) where F<:AbstractFloat
  Thermostat(T, Berendsen!, vars=Berendsen(gamma))
end
```