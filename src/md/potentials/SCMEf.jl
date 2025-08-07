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
SCMEf(; constraints=nothing) = Calculator(SCMEf; EF=SCMEf!, constraints=constraints)

struct _SCMEf_PotVars{T, F<:Float64} <: JMD.PotVars
  lat::Vector{F}
  kwargs::T
end

function SCMEf(bdys::Vector{MyAtoms}) 
  # grab args if user defined them
  kwargs = if isdefined(Main, :SCMEf_USER_ARGS)
    Main.SCMEf_USER_ARGS
  else
    Dict()
  end

  _SCMEf_PotVars(ones(3)*1e3, kwargs)
end

function SCMEf(cell::MyCell)
  lat = reshape(cell.lattice, 9)
  
  # grab args if user defined them
  kwargs = if isdefined(Main, :SCMEf_USER_ARGS)
    Main.SCMEf_USER_ARGS
  else
    Dict()
  end

  _SCMEf_PotVars(lat[[1,5,9]], kwargs)
end

function SCMEf!(F, u, p)
  P = p.potVars
  all(p.PBC) ? pbc = true : pbc = p.PBC

  E, f = py"scmef_get_energy_and_forces"(u, P.lat; pbc=pbc, P.kwargs...)

  F   .= eachrow(f) |> (x -> convert(Vector{Vector{Float64}}, x))

  E
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