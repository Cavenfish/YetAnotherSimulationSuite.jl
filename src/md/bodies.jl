abstract type Particle end

#Atoms in simulation
struct Atom <: Particle
  r::SVector{3,Float64}
  v::SVector{3,Float64}
  m::Float64
  s::Char
end


function getPairs(bdys)

  #distance from vector 
  norm(v) = sqrt(v'v)

  # Get distance matrix 
  d = [[norm(j.r - i.r) for j in bdys] for i in bdys]

  # Get molecules (only works for CO)
  mols = findall.(e -> e < 2, d)

  # Drop duplicates
  mols = unique(mols)
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