# Notes to Self:
#   Making Atom mutable allows mass to swap on fly, but this 
#   might be too niche. Should not cause performance loss and
#   it is nice for the user to have.
#
#   I am making r and v MVectors to restrict array size but still
#   allow mutating the vector. This lets the user translate or add
#   momentum to atoms.

# Atoms in simulation
"""
    Particle{D, F, S}

Mutable struct representing an atom in the simulation.

# Fields
- `r`: Position vector.
- `v`: Velocity vector.
- `m`: Mass.
- `s`: Symbol.
"""
mutable struct Particle{D, F <: AbstractFloat, S <: AbstractString} <: MyAtoms
  r::MVector{D, F}
  v::MVector{D, F}
  m::F
  s::S
end

"""
    Particle(r, v, m, s)

Construct a Particle from position, velocity, mass, and symbol.

# Arguments
- `r`: Position vector.
- `v`: Velocity vector.
- `m`: Mass.
- `s`: Symbol.

# Returns
- Particle object.
"""
function Particle(r::V, v::V, m::Float64, s::String) where V <: Union{Vector, SVector}
  n = length(r)

  rm = MVector{n}(r)
  vm = MVector{n}(v)

  Particle(rm, vm, m, s)
end

"""
    center(bdys::Vector{MyAtoms})

Compute the geometric center of a set of atoms.

# Arguments
- `bdys`: Vector of MyAtoms.

# Returns
- Center position as a vector.
"""
function center(bdys::Vector{MyAtoms})
  o = zeros(3)

  for i in bdys
    o .+= i.r
  end

  o ./ length(bdys)
end

"""
    swapAtoms!(bdys::Vector{MyAtoms}, i, j)

Swap the positions of two atoms in a vector.

# Arguments
- `bdys`: Vector of MyAtoms.
- `i`: Index of first atom.
- `j`: Index of second atom.

# Side Effects
- Modifies the positions in-place.
"""
function swapAtoms!(bdys::Vector{MyAtoms}, i, j)
  a = bdys[i].r |> deepcopy
  b = bdys[j].r |> deepcopy

  bdys[i].r .= b
  bdys[j].r .= a
end

"""
    centerBdys!(bdys::Vector{MyAtoms})

Center a set of atoms at the origin.

# Arguments
- `bdys`: Vector of MyAtoms.

# Side Effects
- Modifies positions in-place.
"""
function centerBdys!(bdys::Vector{MyAtoms})
  com = CoM(bdys)
  for i in bdys
    i.r -= com
  end
end

"""
    swapIso!(bdys::Vector{MyAtoms}, swap, mas)

Swap isotopes (masses) for a set of atoms.

# Arguments
- `bdys`: Vector of MyAtoms.
- `swap`: Indices to swap.
- `mas`: New masses.

# Side Effects
- Modifies masses in-place.
"""
function swapIso!(bdys::Vector{MyAtoms}, swap, mas)
  for i in 1:length(swap)
    j         = swap[i]
    bdys[j].m = mas[i]
  end
end

"""
    getSurfaceMolecules(bdys::Vector{MyAtoms}; α=nothing)

Get the surface molecules from a set of atoms using alpha shapes.

# Arguments
- `bdys`: Vector of MyAtoms.
- `α`: Alpha parameter (optional).

# Returns
- Vector of surface MyAtoms.
"""
function getSurfaceMolecules(bdys::Vector{MyAtoms}; α=nothing)
  mols = getMols(bdys, 1.5)

  pts  = [i.r for i in bdys]

  A    = alphashape(pts; α=α)
  
  i = vcat(A.perimeter...) |> unique
  j = findall.(e -> e in i, mols) |> (x -> findall(e -> !isempty(e), x))

  surf = bdys[vcat(mols[j]...)]

  surf
end

"""
    pickRandomMol(bdys::Vector{MyAtoms}, loc)

Pick a random molecule from the surface or bulk.

# Arguments
- `bdys`: Vector of MyAtoms.
- `loc`: "surf" or "bulk".

# Returns
- Vector of MyAtoms for the selected molecule.
"""
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

"""
    getScaledPos(bdys::Vector{MyAtoms}, lattice)

Get scaled positions for a set of atoms given a lattice.

# Arguments
- `bdys`: Vector of MyAtoms.
- `lattice`: Lattice matrix.

# Returns
- Vector of scaled positions.
"""
function getScaledPos(bdys::Vector{MyAtoms}, lattice)

  T    = inv(lattice)
  
  [T * i.r for i in bdys]
end

"""
    getMIC(bdys::Vector{MyAtoms}, lattice)

Get the minimum image convention positions for a set of atoms.

# Arguments
- `bdys`: Vector of MyAtoms.
- `lattice`: Lattice matrix.

# Returns
- Vector of MyAtoms with minimum image convention applied.
"""
function getMIC(bdys::Vector{MyAtoms}, lattice)
  a,b,c  = eachrow(lattice)
  new    = MyAtoms[]
  s      = repeat([i.s for i in bdys], 27)
  m      = repeat([i.m for i in bdys], 27)
  v      = repeat([i.v for i in bdys], 27)

  # I think it is
  f = [MVector(i*a + j*b + k*c + bdys[q].r)
        for i = -1:1 
          for j = -1:1 
            for k = -1:1 
              for q = 1:length(bdys)]

  for i = 1:length(f)
    push!(new, Particle(f[i], v[i], m[i], s[i]))
  end

  new
end

"""
    wrap!(bdys::Vector{MyAtoms}, lattice)

Wrap all atoms into the primary unit cell.

# Arguments
- `bdys`: Vector of MyAtoms.
- `lattice`: Lattice matrix.

# Side Effects
- Modifies positions in-place.
"""
function wrap!(bdys::Vector{MyAtoms}, lattice)

  f = [floor.(transpose(lattice) \ i.r) for i in bdys]
  
  for i = 1:length(bdys)
    bdys[i].r .-= transpose(lattice) * f[i]
  end

end