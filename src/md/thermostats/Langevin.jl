struct Langevin{F<:AbstractFloat} <:ThermoVars
  gamma::F
  kB::F
end

function Langevin(T::Quantity, gamma::Quantity, calc::MyCalc)
  g = uconvert(calc.time_unit, gamma) |> ustrip

  Thermostat(T, Langevin!, vars=Langevin(g, calc.kB))
end

function Langevin!(a, v, m, Tsim, thermostat)
  N     = length(m)
  gamma = thermostat.vars.gamma
  kB    = thermostat.vars.kB

  # M   = sum(m) / N
  # sig = sqrt(2 * gamma * M * inp.kB * inp.T)
  G   = Normal(0.0, 1)
  eta = [rand(G, 3) for i in 1:N] 

  sigma = @. sqrt(2 * gamma * m * kB * thermostat.T) 
  
  a1 = @. -gamma * v
  a2 = @. sigma / m * eta

  @. a += a1 + a2
end