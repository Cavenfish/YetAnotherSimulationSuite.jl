
# Time dependent MD properties
struct Image{D, F<:AbstractFloat} <: MyImage
  pos::Vector{SVector{D, F}}
  vel::Vector{SVector{D, F}}
  time::F
  temp::F
  energy::F
  forces::Vector{SVector{D, F}}
end

struct Traj{D, F<:AbstractFloat, S<:AbstractString, I<:MyImage} <: MyTraj
  images::Vector{I}
  masses::Vector{F}
  symbols::Vector{S}
  lattice::SMatrix{D,D,F}
end

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

function getBdys(tj::MyTraj, i::Int)
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