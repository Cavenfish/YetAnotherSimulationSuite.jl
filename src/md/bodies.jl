#Atoms in simulation
#Needs mutable to swap masses on the fly
mutable struct Atom <: MyAtoms
  r::Vector{Float64}
  v::Vector{Float64}
  m::Float64
  s::Char
end

function translateBdys!(bdys, v)
  for i in bdys
    i.r .+= v
  end
end

function swapAtoms!(bdys, i, j)
  a = bdys[i].r
  b = bdys[j].r

  bdys[i].r = b
  bdys[j].r = a
end

function centerBdys!(bdys)
  com = CoM(bdys)
  for i in bdys
    i.r -= com
  end
end

function swapIso!(bdys, swap, mas)
  for i in 1:length(swap)
    j         = swap[i]
    bdys[j].m = mas[i]
  end
end

function getSurfaceMolecules(bdys; α=nothing)
  mols = getMols(bdys, 1.5)

  pts  = [i.r for i in bdys]

  A    = alphashape(pts; α=α)
  
  i = vcat(A.perimeter...) |> unique
  j = findall.(e -> e in i, mols) |> (x -> findall(e -> !isempty(e), x))

  surf = bdys[vcat(mols[j]...)]

  surf
end

function pickRandomMol(bdys, loc)

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