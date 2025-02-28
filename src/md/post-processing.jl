
struct Traj <: MyTraj
  t::Vector{Float64}
  E::Vector{Float64}
  T::Vector{Float64}
  F::Vector
  r::Vector
  v::Vector
  m::Vector{Float64}
  s::Vector
end

function processDynamics(solu; dt=fs, step=1)
  t = solu.t[1:step:end] / dt
  N = length(solu.t)

  if typeof(solu.prob.p.energy[1]) == NTuple{5, Float64}
    E = [i[1] for i in solu.prob.p.energy[1:step:end]]
  else
    E = solu.prob.p.energy[1:step:end]
  end

  if typeof(solu.prob.p.forces[1]) == NTuple{5, Float64}
    F = [i[1] for i in solu.prob.p.forces[1:step:end]]
  else
    F = solu.prob.p.forces[1:step:end]
  end

  r = [solu.u[i].x[2] for i in 1:step:N]
  v = [solu.u[i].x[1] for i in 1:step:N]
  m = solu.prob.p.m
  T = [getTemp(m, i, kB, length(m)) for i in v]
  s = [atm.s for atm in solu.prob.p.bdys]

  Traj(t, E, T, F, r, v, m, s)
end

function processDynamics!(tj, solu; dt=fs, step=1)
  t = solu.t[1:step:end] / dt
  N = length(solu.t)

  if typeof(solu.prob.p.energy[1]) == NTuple{5, Float64}
    E = [i[1] for i in solu.prob.p.energy[1:step:end]]
  else
    E = solu.prob.p.energy[1:step:end]
  end

  if typeof(solu.prob.p.forces[1]) == NTuple{5, Float64}
    F = [i[1] for i in solu.prob.p.forces[1:step:end]]
  else
    F = solu.prob.p.forces[1:step:end]
  end

  r = [solu.u[i].x[2] for i in 1:step:N]
  v = [solu.u[i].x[1] for i in 1:step:N]
  m = solu.prob.p.m
  T = [getTemp(m, i, kB, length(m)) for i in v]

  push!(tj.t, t...)
  push!(tj.E, E...)
  push!(tj.T, T...)
  push!(tj.F, F...)
  push!(tj.r, r...)
  push!(tj.v, v...)
end

function processTmpFiles(files; kwargs...)
  
  f    = open(files[1], "r")
  solu = deserialize(f)
  tj   = processDynamics(solu; kwargs...) 
  close(f)

  # Free memory
  @free solu

  for file in files[2:end]
    f    = open(file, "r")
    solu = deserialize(f)
    processDynamics!(tj, solu; kwargs...)
    close(f)

    #Free memory
    @free solu
  end

  return tj
end