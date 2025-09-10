struct Berendsen{F<:AbstractFloat} <:ThermoVars
  gamma::F
end

function Berendsen(T::Quantity, gamma::Quantity, calc::MyCalc)
  g = uconvert(calc.time_unit, gamma) |> ustrip

  Thermostat(T, Berendsen!, vars=Berendsen(g))
end

function Berendsen!(a, v, m, Tsim, thermostat)
  gamma = thermostat.vars.gamma

  if Tsim == 0.0
    a .+= gamma .* v
  else
    a .+= gamma * (thermostat.T / Tsim - 1) .* v
  end
end