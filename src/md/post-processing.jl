
struct Traj
  t::Vector{Float64}
  E::Vector{Float64}
  T::Vector{Float64}
  F::Vector
  sys::Vector
end

function processDynamics(solu)
  t = solu.t

  if typeof(solu.prob.p.energy[1]) == Dict{Any, Any}
    E = [i["total"] for i in solu.prob.p.energy]
  else
    E = [i for i in solu.prob.p.energy]
  end

  if typeof(solu.prob.p.forces[1]) == Dict{Any, Any}
    F = [i["total"] for i in solu.prob.p.forces]
  else
    F = [i for i in solu.prob.p.forces]
  end

  v,m = getVelMas(solu)
  T   = [getTemp(m, i, kB, length(m)) for i in v]

  sys  = []
  bdys = solu.prob.p.bdys
  for i in 1:length(t)
    frame = Atom[]
    for j in 1:length(bdys)
      r = solu.u[i].x[2][j]
      v = solu.u[i].x[1][j]
      m = bdys[j].m
      s = bdys[j].s
      push!(frame, Atom(r,v,m,s))
    end
    push!(sys, frame)
  end

  return Traj(t, E, T, F, sys)
end
