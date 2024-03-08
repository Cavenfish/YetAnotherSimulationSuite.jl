
#Atoms in simulation
mutable struct Atom
  r::Vector{Float64}
  v::Vector{Float64}
  m::Float64
  s::Char
end

struct Molecule{N}
  r::SMatrix{N,3,Float64}
  v::SMatrix{N,3,Float64}
  m::SVector{N,Float64}
  name::String
end

function getMols(bdys, rmin; D=3)
  r   = [i.r for i in bdys]

  pts = hcat(r...)

  ret = dbscan(pts[1:D, :], rmin)

  [i.core_indices for i in ret.clusters]
end

function getPairs(bdys)

  # Get mols and N
  mols = if length(bdys) < 3
    getMols(bdys, 1.5, D=length(bdys)-1) 
  else
    getMols(bdys, 1.5)
  end
  N    = size(mols)[1]

  # Make all pairs
  pars = Pair[]
  for i in 1:N
    for j in i+1:N
      push!(pars, Pair(mols[i],mols[j]))
    end
  end

  return pars, mols
end