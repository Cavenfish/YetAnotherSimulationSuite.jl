"""
    processDynamics(solu::SciMLBase.ODESolution; dt=fs, step=1)

Convert an ODE solution to a Traj object.

# Arguments
- `solu`: ODE solution object.
- `dt`: Time step (default: fs).
- `step`: Step interval (default: 1).

# Returns
- Traj object.
"""
function processDynamics(solu::SciMLBase.ODESolution; dt=fs, step=1)
  N  = length(solu.t)
  
  Traj(
    [getImage(solu, i, dt) for i = 1:step:N],
    solu.prob.p.m,
    solu.prob.p.s,
    solu.prob.p.ensemble.lattice
  )
end

"""
    push!(tj::MyTraj, solu::SciMLBase.ODESolution; dt=fs, step=1)

Append images from an ODE solution to a trajectory.

# Arguments
- `tj`: Traj object.
- `solu`: ODE solution object.
- `dt`: Time step (default: fs).
- `step`: Step interval (default: 1).

# Side Effects
- Modifies `tj` in-place.
"""
function Base.push!(tj::MyTraj, solu::SciMLBase.ODESolution; dt=fs, step=1)
  N = length(solu.t)
  
  for i = 1:step:N
    push!(tj.images, getImage(solu, i, dt))
  end

end

"""
    processTmpFiles(files; kwargs...)

Process a list of temporary files containing ODE solutions into a single trajectory.

# Arguments
- `files`: List of file paths.
- `kwargs`: Additional keyword arguments.

# Returns
- Traj object.
"""
function processTmpFiles(files; kwargs...)
  
  # Process first tmp file into traj obj
  tj = open(files[1], "r") do io
    solu = deserialize(io)
    processDynamics(solu; kwargs...)
  end

  # Free memory. Needed?
  @free solu

  for file in files[2:end]

    # Append traj with next tmp file data
    open(file, "r") do io
      solu = deserialize(io)
      push!(tj, solu; kwargs...)
    end

    #Free memory. Needed?
    @free solu
  end

  # Clean-up
  for f in files
    rm(f)
  end

  tj
end