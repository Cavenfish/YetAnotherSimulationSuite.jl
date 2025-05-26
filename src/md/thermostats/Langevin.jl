struct Langevin{F<:AbstractFloat} <:ThermoVars
  gamma::F
end

function Langevin(T::F, gamma::F) where F<:AbstractFloat
  Thermostat(T, Langevin!, vars=Langevin(gamma))
end

function Langevin!(a, v, m, Tsim, thermostat)
  N     = length(m)
  gamma = thermostat.vars.gamma

  # M   = sum(m) / N
  # sig = sqrt(2 * gamma * M * inp.kB * inp.T)
  G   = Normal(0.0, 1)
  eta = [rand(G, 3) for i in 1:N] 

  sigma = @. sqrt(2 * gamma * m * kB * thermostat.T) 
  
  a1 = @. -gamma * v
  a2 = @. sigma / m * eta

  @. a += a1 + a2
end