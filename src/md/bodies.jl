# Notes to Self:
#   Making Atom mutable allows mass to swap on fly, but this 
#   might be too niche. Most users likely don't need that, so I
#   am scrubbing it.
#
#   I am making r and v MVectors to restrict array size but still
#   allow mutating the vector. This lets the user translate or add
#   momentum to atoms.

# Atoms in simulation
struct Particle{D} <: MyAtoms
  r::MVector{D, Float64}
  v::MVector{D, Float64}
  m::Float64
  s::String
end

function translateBdys!(bdys::Vector{MyAtoms}, v)
  for i in bdys
    i.r .+= v
  end
end

function swapAtoms!(bdys::Vector{MyAtoms}, i, j)
  a = bdys[i].r
  b = bdys[j].r

  bdys[i].r = b
  bdys[j].r = a
end

function centerBdys!(bdys::Vector{MyAtoms})
  com = CoM(bdys)
  for i in bdys
    i.r -= com
  end
end

function swapIso!(bdys::Vector{MyAtoms}, swap, mas)
  for i in 1:length(swap)
    j         = swap[i]
    bdys[j].m = mas[i]
  end
end

function getSurfaceMolecules(bdys::Vector{MyAtoms}; α=nothing)
  mols = getMols(bdys, 1.5)

  pts  = [i.r for i in bdys]

  A    = alphashape(pts; α=α)
  
  i = vcat(A.perimeter...) |> unique
  j = findall.(e -> e in i, mols) |> (x -> findall(e -> !isempty(e), x))

  surf = bdys[vcat(mols[j]...)]

  surf
end

function pickRandomMol(bdys::Vector{MyAtoms}, loc)

  surf = getSurfaceMolecules(bdys)
  bulk = [i for i in bdys if !(i in surf)]

  if loc == "surf"

    mols = getMols(bdys, 1.5)
    return surf[rand(mols)]

  elseif loc == "bulk"

    mols = getMols(bdys, 1.5)
    return bulk[rand(mols)]

  end
  
end

function getScaledPos(bdys::Vector{MyAtoms}, lattice)

  T    = inv(lattice)
  
  [T * i.r for i in bdys]
end

function getMIC(bdys::Vector{MyAtoms}, lattice)
  a,b,c  = eachrow(lattice)
  new    = MyAtoms[]
  s      = repeat([i.s for i in bdys], 27)
  m      = repeat([i.m for i in bdys], 27)
  v      = repeat([i.v for i in bdys], 27)

  # I think it is
  f = [i*a + j*b + k*c + bdys[q].r
        for i = -1:1 
          for j = -1:1 
            for k = -1:1 
              for q = 1:length(bdys)]

  for i = 1:length(f)
    push!(new, Atom(f[i], v[i], m[i], s[i]))
  end

  new
end

function wrap!(bdys::Vector{MyAtoms}, lattice)

  f = [floor.(transpose(lattice) \ i.r) for i in bdys]
  
  for i = 1:length(bdys)
    bdys[i].r .-= transpose(lattice) * f[i]
  end

end