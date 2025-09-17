#TODO:
#   -Make function to clean duplicates

"""
    Cell{D, B, I, F, S}

Structure representing a simulation cell.

# Fields
- `lattice`: Lattice matrix.
- `scaled_pos`: Scaled positions.
- `velocity`: Velocities.
- `masses`: Masses.
- `symbols`: Atomic symbols.
- `mask`: Mask vector.
- `PBC`: Periodic boundary conditions.
- `NC`: Neighbor counts.
"""
struct Cell{D, B, I<:Int, F<:AbstractFloat, S<:AbstractString} <: MyCell
  lattice::MMatrix{D,D,F}
  scaled_pos::Vector{MVector{D,F}}
  velocity::Vector{MVector{D,F}}
  masses::Vector{F}
  symbols::Vector{S}
  mask::Vector{B}
  PBC::Vector{B}
  NC::Vector{I}
end

"""
    Cell(lat, spos, vel, mas, sym, mask, PBC, NC)

Construct a Cell object from lattice, positions, velocities, etc.

# Arguments
- `lat`: Lattice matrix.
- `spos`: Scaled positions.
- `vel`: Velocities.
- `mas`: Masses.
- `sym`: Symbols.
- `mask`: Mask vector.
- `PBC`: Periodic boundary conditions.
- `NC`: Neighbor counts.

# Returns
- Cell object.
"""
function Cell(
  lat::Matrix, spos::AbstractArray, vel::AbstractArray, 
  mas::Vector{F}, sym::Vector{S}, mask::Vector{Bool},
  PBC::Vector{Bool}, NC::Vector{Int}
) where {F<:AbstractFloat, S<:AbstractString}

  n = length(spos[1])

  Cell(
    MMatrix{n,n}(lat),
    MVector{n}.(spos),
    MVector{n}.(vel),
    mas, sym, mask, PBC, NC
  )
end

"""
    length(cell::MyCell)

Get number of atoms in cell

# Arguments
- `cell`: MyCell object

# Returns
- Number of atoms in cell
"""
Base.length(cell::MyCell) = length(cell.masses)

"""
    show(io::IO, cell::MyCell)

Display cell struct information

# Arguments
- `io`: IO 
- `cell`: MyCell object
"""
function Base.show(io::IO, cell::MyCell)
  s = unique(cell.symbols)
  N = length(cell.masses)
  println(io, "$(N) Atoms")
  
  for i in s
    n = filter(e -> e == i, cell.symbols) |> length
    println(io, "$(n) $(i)")
  end

  println(io, "lattice:\n$(cell.lattice)")
end

"""
    getindex(cell::MyCell, inds)

Get cell atoms at `inds` indicies

# Arguments
- `cell`: MyCell object
- `inds`: Indices to get

# Returns
- Cell object
"""
function Base.getindex(cell::MyCell, inds)
  Cell(
    copy(cell.lattice),
    cell.scaled_pos[inds],
    cell.velocity[inds],
    cell.masses[inds],
    cell.symbols[inds],
    cell.mask[inds],
    copy(cell.PBC),
    copy(cell.NC)
  )
end

"""
    deleteat!(cell::MyCell, iter)

Remove atoms at given indices from a cell.

# Arguments
- `cell`: MyCell object.
- `iter`: Indices to remove.

# Side Effects
- Modifies the cell in-place.
"""
function Base.deleteat!(cell::MyCell, iter)
  deleteat!(cell.scaled_pos, iter)
  deleteat!(cell.velocity, iter)
  deleteat!(cell.masses, iter)
  deleteat!(cell.symbols, iter)
  deleteat!(cell.mask, iter)
end

"""
    reorder!(cell::MyCell, order::Vector{Int})

Reorder the atoms in a cell according to a given order.

# Arguments
- `cell`: MyCell object.
- `order`: New order of indices.

# Side Effects
- Modifies the cell in-place.
"""
function reorder!(cell::MyCell, order::Vector{Int})
  cell.scaled_pos .= cell.scaled_pos[order]
  cell.velocity   .= cell.velocity[order]
  cell.masses     .= cell.masses[order]
  cell.symbols    .= cell.symbols[order]
  cell.mask       .= cell.mask[order]
end

"""
    makeCell(bdys, lattice; mask, PBC, NC)

Construct a Cell from a vector of MyAtoms and a lattice.

# Arguments
- `bdys`: Vector of MyAtoms.
- `lattice`: Lattice matrix.
- `mask`: Mask vector (optional).
- `PBC`: Periodic boundary conditions (optional).
- `NC`: Neighbor counts (optional).

# Returns
- Cell object.
"""
function makeCell(bdys::Vector{MyAtoms}, lattice::AbstractMatrix; 
  mask=repeat([false], length(bdys)), PBC=repeat([true], 3), NC=[1,1,1])

  Cell(
    lattice,
    getScaledPos(bdys, lattice),
    [i.v for i in bdys],
    [i.m for i in bdys],
    [i.s for i in bdys],
    mask, PBC, NC
  )
end

"""
    makeBdys(cell::MyCell)::Vector{MyAtoms}

Convert a cell to a vector of MyAtoms.

# Arguments
- `cell`: MyCell object.

# Returns
- Vector of MyAtoms.
"""
function makeBdys(cell::MyCell)::Vector{MyAtoms}
  D    = length(cell.velocity[1])
  pos  = MVector{D}.(getPos(cell))
  vel  = [MVector{D}(i) for i in cell.velocity]

  [Particle(pos[i], vel[i], cell.masses[i], cell.symbols[i]) for i = 1:length(pos)]
end

"""
    getScaledPos(x0, lattice)

Get scaled positions from Cartesian coordinates and lattice.

# Arguments
- `x0`: Cartesian coordinates.
- `lattice`: Lattice matrix.

# Returns
- Vector of scaled positions.
"""
function getScaledPos(x0, lattice)
  T = inv(lattice)
  
  [T * x0[i:i+2] for i = 1:3:length(x0)-1]
end

"""
    getScaledPos!(cell, x0)

Update scaled positions in a cell from Cartesian coordinates.

# Arguments
- `cell`: MyCell object.
- `x0`: Cartesian coordinates.

# Side Effects
- Modifies the cell in-place.
"""
function getScaledPos!(cell, x0)
  T = inv(cell.lattice)

  for i = 1:3:length(x0)-1
    j::Int = (i+2)/3

    cell.scaled_pos[j] .= T * x0[i:i+2]
  end

end

"""
    getPos(cell)

Get Cartesian positions from scaled positions in a cell.

# Arguments
- `cell`: MyCell object.

# Returns
- Vector of Cartesian positions.
"""
function getPos(cell)
  [cell.lattice * i for i in cell.scaled_pos]
end

"""
    getVolume(cell)

Compute the volume of a cell.

# Arguments
- `cell`: MyCell object.

# Returns
- Volume (Float64).
"""
function getVolume(cell)
  a,b,c  = eachrow(cell.lattice)

  cross(b,c) |> (x -> dot(a,x))
end 

"""
    wrap!(cell)

Wrap all atoms in a cell into the primary unit cell.

# Arguments
- `cell`: MyCell object.

# Side Effects
- Modifies the cell in-place.
"""
function wrap!(cell)
  r = getPos(cell)
  
  f = [floor.(transpose(cell.lattice) \ i) for i in r]

  for i = 1:length(r)
    r[i] .-= transpose(cell.lattice) * f[i]
  end

  T = inv(cell.lattice)

  for i = 1:length(r)
    cell.scaled_pos[i] .= T * r[i]
  end

end

"""
    center!(cell)

Center the atoms in a cell.

# Arguments
- `cell`: MyCell object.

# Side Effects
- Modifies the cell in-place.
"""
function center!(cell)
  r  = getPos(cell)
  r0 = sum(r) ./ length(r)
  μ  = sum(0.5 * cell.lattice, dims=1) |> (x -> reshape(x, 3))
  x  = μ - r0
  T  = inv(cell.lattice)
  
  # Translate positions and scale positions
  for i in r
    i .+= x
    i .= T * i
  end

  # Update scaled positions
  for i in 1:length(r)
    cell.scaled_pos[i] .= r[i]
  end 

end

"""
    repeat(cell::MyCell, count::Integer)

Repeat a cell multiple times.

# Arguments
- `cell`: MyCell object.
- `count`: Number of repetitions.

# Returns
- New repeated Cell object.
"""
function Base.repeat(cell::MyCell, count::Integer)
  Cell(
    deepcopy(cell.lattice),
    myRepeat(cell.scaled_pos, count, cell.mask),
    myRepeat(cell.velocity, count, cell.mask),
    myRepeat(cell.masses, count, cell.mask),
    myRepeat(cell.symbols, count, cell.mask),
    myRepeat(cell.mask, count, cell.mask),
    deepcopy(cell.PBC), 
    deepcopy(cell.NC)
  )
end

"""
    replicate!(super, cell, N)

Replicate a cell into a supercell.

# Arguments
- `super`: Supercell object to fill.
- `cell`: Cell to replicate.
- `N`: Replication factors along each axis.

# Side Effects
- Modifies the supercell in-place.
"""
function replicate!(super, cell, N)

  a,b,c  = eachrow(cell.lattice)
  m      = length(cell.symbols)
  M      = length(super.scaled_pos)
  n      = prod(N) * m
  l      = n - (sum(cell.mask) * (prod(N) - 1))

  # End if super is too small 
  # TODO: make this increase super size
  if M < l 
    @error "super is too small"
    return
  end

  # Trim super (if needed)
  if M > l
    trim!(super, l+1:M)
  end

  x = 1
  for i = 0:N[1]-1
    for j = 0:N[2]-1
      for k = 0:N[3]-1
        for q = 1:m
          # Skip replication if atom is masked
          cell.mask[q] && x > m && continue

          # Get replicated pos
          r   = cell.lattice * cell.scaled_pos[q]
          r .+= (i*a) + (j*b) + (k*c)

          # inplace update super
          super.scaled_pos[x] .= inv(super.lattice) * r
          super.velocity[x]   .= cell.velocity[q]
          super.masses[x]      = cell.masses[q]
          super.symbols[x]     = cell.symbols[q]
          super.mask[x]        = cell.mask[q]
          
          # increment x
          x += 1
        end
      end
    end
  end

end

"""
    makeSuperCell(cell, T)

Create a supercell from a cell and transformation matrix.

# Arguments
- `cell`: MyCell object.
- `T`: Transformation matrix.

# Returns
- New supercell object.
"""
function makeSuperCell(cell, T)

  N       = diag(T)
  lattice = T * cell.lattice

  # Clean up machine precision noise
  for i = 1:9
    if abs(lattice[i]) < 1e-8
      lattice[i] = 0.0
    end 
  end

  # allocate super cell
  super = repeat(cell, prod(N))

  # inplace update lattice
  super.lattice .= lattice

  # inplace update fields
  replicate!(super, cell, N)

  # inplace wrap atoms outside PBC
  wrap!(super)

  super
end

"""
    makeSuperCell!(super, cell, T)

In-place creation of a supercell from a cell and transformation matrix.

# Arguments
- `super`: Supercell object to fill.
- `cell`: Cell to replicate.
- `T`: Transformation matrix.

# Side Effects
- Modifies the supercell in-place.
"""
function makeSuperCell!(super, cell, T)

  super.lattice .= T * cell.lattice

  # Clean up machine precision noise
  for i = 1:9
    if abs(super.lattice[i]) < 1e-8
      super.lattice[i] = 0.0
    end 
  end
    
  replicate!(super, cell, diag(T))
  wrap!(super)

end