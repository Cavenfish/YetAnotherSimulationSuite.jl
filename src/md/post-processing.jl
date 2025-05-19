
function processDynamics(solu::SciMLBase.ODESolution; dt=fs, step=1)
  N  = length(solu.t)
  s  = [i.s for i in solu.prob.p.bdys] # WILL CHANGE
  tj = Traj(
    getImage(solu, 1, dt),
    solu.prob.p.m,
    s,
    solu.prob.p.lattice
  )
  
  for i = 2:N
    push!(tj.images, getImage(solu, i, dt))
  end

  tj
end

function Base.push!(tj::MyTraj, solu::SciMLBase.ODESolution; dt=fs, step=1)
  N = length(solu.t)
  
  for i = 1:N
    push!(tj.images, getImage(solu, i, dt))
  end
  
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