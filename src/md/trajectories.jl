
# Time dependent MD properties
struct Image{D, F<:AbstractFloat} <: MyImage
  pos::Vector{SVector{D, F}}
  vel::Vector{SVector{D, F}}
  time::F
  temp::F
  energy::F
  forces::Vector{SVector{D, F}}
end

struct Traj{D, F<:AbstractFloat, S<:AbstractString, Im<:MyImage} <: MyTraj
  images::Vector{Im}
  masses::Vector{F}
  symbols::Vector{S}
  lattice::SMatrix{D,D,F}
end

"""
    length(traj::MyTraj)

Get the number of images in a trajectory.

# Arguments
- `traj`: Traj object.

# Returns
- Number of images.
"""
Base.length(traj::MyTraj) = length(traj.images)

"""
    show(io::IO, traj::MyTraj)

Custom display for MyTraj objects.

# Arguments
- `io`: IO stream.
- `traj`: MyTraj object.
"""
function Base.show(io::IO, traj::MyTraj)
  N = length(traj)
  println(io, "$(N) Images")
end

"""
    Traj(imgs, mas, sym, lat)

Construct a Traj object from images, masses, symbols, and lattice.

# Arguments
- `imgs`: Vector of Image objects.
- `mas`: Vector of masses.
- `sym`: Vector of symbols.
- `lat`: Lattice matrix.

# Returns
- Traj object.
"""
function Traj(
  imgs::T, mas::Vector{Float64}, 
  sym::Vector{String}, lat::Union{Matrix, MMatrix}
) where T

  n,_     = size(lat)
  lattice = SMatrix{n,n}(lat)

  Traj(imgs, mas, sym, lattice)
end

"""
    getImage(solu::SciMLBase.ODESolution, i::Int, dt::Float64)

Extract an Image from an ODE solution at a given index.

# Arguments
- `solu`: ODE solution object.
- `i`: Index of the time step.
- `dt`: Time step size.

# Returns
- Image object.
"""
function getImage(solu::SciMLBase.ODESolution, i::Int, dt::Float64)
  r = solu.u[i].x[2]
  v = solu.u[i].x[1]
  t = solu.t[i] / dt
  T = solu.prob.p.temp[i]
  E = solu.prob.p.energy[i]
  n = length(r[1])
  F = SVector{n}.(solu.prob.p.forces[i])

  Image(r, v, t, T, E, F)
end

"""
    makeBdys(tj::MyTraj, i::Int)

Construct a vector of MyAtoms from a trajectory at a given image index.

# Arguments
- `tj`: Traj object.
- `i`: Image index.

# Returns
- Vector of MyAtoms.
"""
function makeBdys(tj::MyTraj, i::Int)
  bdys = MyAtoms[]

  for j = 1:length(tj.masses)
    Particle(
      tj.images[i].pos[j],
      tj.images[i].vel[j],
      tj.masses[j],
      tj.symbols[j]
    ) |> (x -> push!(bdys, x))
  end
  
  bdys
end

"""
    getVelMas(tj::MyTraj)

Get all velocities and masses from a trajectory.

# Arguments
- `tj`: Traj object.

# Returns
- Tuple: (vector of velocities, vector of masses)
"""
function getVelMas(tj::MyTraj)
  vel = [i.vel for i in tj.images]
  
  vel, tj.masses
end

function getLastFrame!(bdys::Vector{MyAtoms}, tj::MyTraj)
  n = length(tj)

  for i = 1:length(bdys)
    bdys[i].r .= tj.images[n].pos[i]
    bdys[i].v .= tj.images[n].vel[i]
  end
end