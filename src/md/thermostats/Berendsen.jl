struct Berendsen{F<:AbstractFloat} <:ThermoVars
  gamma::F
end

function Berendsen(T::F, gamma::F) where F<:AbstractFloat
  Thermostat(T, Berendsen!, vars=Berendsen(gamma))
end

function Berendsen!(a, v, m, Tsim, thermostat)
  gamma = thermostat.vars.gamma

  if Tsim == 0.0
    a .+= gamma .* v
  else
    a .+= gamma * (thermostat.T / Tsim - 1) .* v
  end
end