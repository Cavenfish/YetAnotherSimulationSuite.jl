"""
SCME/f potential 

This calls an ASE calculator of the SCME/f potential (not great)
I still want to make a better wrapper. Idealy, one that calls
the c++ routines directly (like MBX).

SCME Paper:
https://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp52097h

SCME/f Paper:
https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00598
"""

struct _SCMEf_PotVars <: JMD.PotVars
  lat
end

function SCMEf(cell::JMD.MyCell)
  lat = reshape(cell.lattice, 9)

  _SCMEf_PotVars(lat[[1,5,9]])
end

function SCMEf(dv, v, u, p, t)
  all(p.PBC) ? pbc = true : pbc = p.PBC

  E, f = py"scmef_get_energy_and_forces"(u, p.potVars.lat, NC=p.NC, pbc=pbc)

  F    = eachrow(f) |> (x -> convert(Vector{Vector{Float64}}, x))

  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end

function SCMEf(F, G, y0, p)

  # initialize things
  E = 0.0
  u = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
  end

  all(p.PBC) ? pbc = true : pbc = p.PBC

  if G != nothing
    e, f = py"scmef_get_energy_and_forces"(u, p.potVars.lat, NC=p.NC, pbc=pbc)
    E    = e
    G   .= transpose(-f) |> (x -> reshape(x, length(x)))
  else
    E    = py"scmef_get_energy"(u, p.potVars.lat, NC=p.NC, pbc=pbc)
  end

  if F != nothing
    return E
  end

end

function SCMEf(F, G, cell::JMD.MyCell, lat)
  tmp          = deepcopy(cell)
  tmp.lattice .= reshape(lat, (3,3))

  # Force all bond lengths to be 0.975 to prevent bias
  pos  = getPos(cell)
  for i = 1:3:length(pos)
    ro, rh1, rh2 = pos[[i, i+1, i+2]]
    
    rvec  = ro - rh1
    r     = norm(rvec)
    x     = 0.975 - r
    rh1 .-= x * (rvec / r)

    rvec  = ro - rh2
    r     = norm(rvec)
    x     = 0.975 - r
    rh2 .-= x * (rvec / r)

    #Maybe I can condense this somehow
  end

  # Orthogonal stress because SCME/f only does orthogonal
  if G != nothing
    G .= getNumericalStressOrthogonal(SCMEf, tmp) |> (x -> reshape(x, 9))
  end

  if F != nothing
    return getPotEnergy(SCMEf, tmp)
  end
end