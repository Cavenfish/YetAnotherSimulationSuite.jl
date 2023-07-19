
struct Traj
  t::Vector{Float64}
  E::Vector{Float64}
  T::Vector{Float64}
  F::Vector
  r::Vector
  v::Vector
  m::Vector{Float64}
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

  return Traj(t, E, T, F, r, v, m)
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

  # Debating if I want to add this
  # f    = open("$(file[1]).v", "w")
  # v    = [solu.u[i].x[1] for i in 1:length(solu.t)]
  # serialize(f, v)
  # close(f)

  # Free memory
  @free solu

  for file in files[2:end]
    f    = open(file, "r")
    solu = deserialize(f)
    processDynamics!(tj, solu; kwargs...)
    close(f)

    # Debating if I want to add this
    # f    = open("$file.v", "w")
    # v    = [solu.u[i].x[1] for i in 1:length(solu.t)]
    # serialize(f, v)
    # close(f)

    #Free memory
    @free solu
    # @free v
  end

  return tj
end

function trackEnergyDissipation(traj, pot, mol)

  m = traj.m
  N = length(traj.t)
  n = length(m)

  molVib, molRot, molTra = [],[],[]
  avgVib, avgRot, avgTra = [],[],[]

  others = [[i,i+1] for i in 1:2:n if !(i in mol)]

  for i in 1:N
    r = traj.r[i]
    v = traj.v[i]

    mvib = getCOVibEnergy(r[mol], v[mol], m[mol]; pot=pot)
    mrot = getRotEnergy(  r[mol], v[mol], m[mol])
    mtra = getTransEnergy(r[mol], v[mol], m[mol])

    push!(molVib, mvib)
    push!(molRot, mrot)
    push!(molTra, mtra)

    avib,arot,atra = 0.0,0.0,0.0
    for j in others
      avib += getCOVibEnergy(r[j], v[j], m[j]; pot=pot)
      arot += getRotEnergy(  r[j], v[j], m[j])
      atra += getTransEnergy(r[j], v[j], m[j])
    end


    push!(avgVib, avib/n)
    push!(avgRot, arot/n)
    push!(avgTra, atra/n)
  end

  df = DataFrame(  time=traj.t,
                   temp=traj.T,
                 molVib=molVib,
                 molRot=molRot,
                 molTra=molTra,
                 avgVib=avgVib,
                 avgRot=avgRot,
                 avgTra=avgTra)

  return df
end
