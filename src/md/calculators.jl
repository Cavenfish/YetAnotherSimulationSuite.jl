
struct Calculator{B, E, F, EF} <: MyCalc
  b::B
  e::E
  f!::F
  ef!::EF
end

function Calculator(b; E=nothing, F=nothing, EF=nothing)
  if E == nothing && EF == nothing
    throw()
  end

  if F == nothing && EF == nothing
    throw()
  end

  Calculator(b, E, F, EF)
end

function fg!(F, G, x, p, calc::MyCalc)
  # initialize things
  u      = [x[i:i+2] for i = 1:3:length(x)]
  forces = [zeros(3) for i = 1:3:length(x)]

  if F != nothing && G != nothing
    
    E = if calc.ef! != nothing
      calc.ef!(forces, u, p)
    else
      calc.f!(forces, u, p)
      calc.e(u, p)
    end

    tmp = [j for i in forces for j in i]
    for i in 1:length(G)
      G[i] = -tmp[i]
    end

    return E
  end

  if F != nothing 
    if calc.e != nothing
      return calc.e(u, p)
    else
      return calc.ef!(forces, u, p)
    end
  end

  if G != nothing

    if calc.f! != nothing
      calc.f!(forces, u, p)
    else
      calc.ef!(forces, u, p)
    end

    tmp = [j for i in forces for j in i]
    for i in 1:length(G)
      G[i] = -tmp[i]
    end
  end

end

function dyn!(dv, v, u, p, t, calc)
  # initialize things
  F = [@MVector zeros(3) for i = 1:length(u)]

  # Calculate energy and forces
  E = if calc.ef! != nothing
    calc.ef!(F, u, p)
  else
    calc.f!(F, u, p)
    calc.e(u, p)
  end

  dv .= F ./ p.m
  T   = getTemp(p.m, v)

  if isa(p.ensemble, NVT)
    thermostat = p.ensemble.thermostat
    thermostat.act!(dv, v, p.m, T, thermostat)
  end

  push!(p.temp,   T)
  push!(p.energy, E)
  push!(p.forces, F)
end