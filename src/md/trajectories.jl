

struct Image{D} <: MyImage
  pos::Vector{SVector{D, Float64}}
  vel::Vector{SVector{D, Float64}}
  time::Float64
  temp::Float64
  energy::Float64
  forces::Vector{SVector{D, Float64}}
end

struct Traj <: MyTraj
  images::Vector{MyImage}
  masses::Vector{Float64}
  symbols::Vector{Float64}
  lattice::Matrix{Float64}
end

function getImage(solu::SciMLBase.ODESolution, i::Int, dt::Float64)
  r = solu.u[i].x[2]
  v = solu.u[i].x[1]
  t = solu.t[i] / dt
  T = solu.prob.p.temp[i]
  E = solu.prob.p.energy[i]
  F = solu.prob.p.forces[i]

  Image(r, v, t, T, E, F)
end

function getLastImage(solu)
  n = length(solu.prob.p.bdys)
  
  new = Atom[]
  for i in 1:n
    r = solu.u[end].x[2][i]
    v = solu.u[end].x[1][i]
    m = solu.prob.p.bdys[i].m
    s = solu.prob.p.bdys[i].s
    push!(new, Particle(r, v, m, s))
  end
  return new
end

function getLastImage!(cell::MyCell, solu)
  x0 = [j for i in solu.u[end].x[2] for j in i]

  cell.velocity   .= solu.u[end].x[1]
  cell.scaled_pos .= getScaledPos(x0, cell.lattice)

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