# TODO:
#   - Generalize the dimensionality of optimizations

"""
    optVars

Structure holding optimization variables for geometry optimization.

# Fields
- `potVars`: Potential variables.
- `mols`: Molecule indices.
- `m`: Masses.
- `PBC`: Periodic boundary conditions.
- `NC`: Neighbor counts.
- `lattice`: Lattice matrix.
"""
struct optVars{B, I<:Int, PV<:PotVars, F<:AbstractFloat, AM<:AbstractMatrix}
  potVars::PV
  mols::Vector{Vector{I}}
  m::Vector{F}
  PBC::Vector{B}
  NC::Vector{I}
  lattice::AM
end

"""
    prepX0(bdys::Vector{MyAtoms})

Prepare the initial coordinate vector from a set of atoms.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.

# Returns
- Flat vector of coordinates.
"""
function prepX0(bdys::Vector{MyAtoms})
  r  = [i.r for i in bdys]
  
  [j for i in r for j in i]
end

"""
    prepX0(cell::MyCell)

Prepare the initial coordinate vector from a cell.

# Arguments
- `cell`: `MyCell` object.

# Returns
- Flat vector of coordinates.
"""
function prepX0(cell::MyCell)
  r  = getPos(cell)
  
  [j for i in r for j in i]
end

"""
    prep4pot(builder, bdys::Vector{MyAtoms})

Prepare variables for potential energy calculation from atoms.

# Arguments
- `builder`: Potential builder function.
- `bdys`: Vector of `MyAtoms` objects.

# Returns
- Tuple: (coordinate vector, `optVars` object).
"""
function prep4pot(builder, bdys::Vector{MyAtoms})
  m       = [i.m for i in bdys]
  x0      = prepX0(bdys)
  potVars = builder(bdys)
  mols    = getMols(bdys, 1.2)
  NC      = [0,0,0]
  PBC     = repeat([false], 3)
  lattice = MMatrix{3,3}(zeros(3,3))
  vars    = optVars(potVars, mols, m, PBC, NC, lattice)
  
  x0, vars
end

"""
    prep4pot(builder, cell::MyCell)

Prepare variables for potential energy calculation from a cell.

# Arguments
- `builder`: Potential builder function.
- `cell`: `MyCell` object.

# Returns
- Tuple: (coordinate vector, `optVars` object).
"""
function prep4pot(builder, cell::MyCell)
  bdys    = makeBdys(cell)
  x0      = prepX0(cell)
  potVars = builder(cell)
  mols    = getMols(cell, 1.2)
  vars    = optVars(potVars, mols, cell.masses, 
                    cell.PBC, cell.NC, cell.lattice)
  
  x0, vars
end

"""
    getNewBdys(bdys, res)

Construct new atom objects from optimization results.

# Arguments
- `bdys`: Original vector of `MyAtoms`.
- `res`: Optimization result.

# Returns
- Vector of new `MyAtoms` objects.
"""
function getNewBdys(bdys, res)
  N   = length(bdys)
  opt = res.minimizer
  new = MyAtoms[]

  for i in 1:3:N*3
    j::UInt16 = (i+2) / 3

    r = MVector{3}(opt[i:i+2])
    v = bdys[j].v
    m = bdys[j].m
    s = bdys[j].s
    push!(new, Particle(r,v,m,s))
  end

  new
end

"""
    opt(calc::MyCalc, algo, bdys::Vector{MyAtoms}; kwargs...)

Optimize the geometry of a set of atoms.

# Arguments
- `calc`: Calculator object.
- `algo`: Optimization algorithm.
- `bdys`: Vector of `MyAtoms` objects.
- `kwargs`: Additional options.

# Returns
- Optimized vector of `MyAtoms`.
"""
function opt(calc::MyCalc, algo, bdys::Vector{MyAtoms}; kwargs...)
  x0, vars = prep4pot(calc.b, bdys)
  optFunc  = NLSolversBase.only_fg!((F,G,x) -> fg!(F,G,x, vars, calc))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  optBdys  = getNewBdys(bdys, res)

  optBdys
end

"""
    opt(calc::MyCalc, algo, cell::MyCell; kwargs...)

Optimize the geometry of a cell.

# Arguments
- `calc`: Calculator object.
- `algo`: Optimization algorithm.
- `cell`: `MyCell` object.
- `kwargs`: Additional options.

# Returns
- Optimized `MyCell` object.
"""
function opt(calc::MyCalc, algo, cell::MyCell; kwargs...)
  x0, vars = prep4pot(calc.b, cell)
  optFunc  = NLSolversBase.only_fg!((F,G,x) -> fg!(F,G,x, vars, calc))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  spos     = getScaledPos(res.minimizer, cell.lattice)
  spos     = [MVector(i) for i in spos]

  # The new cell returned will have fields that point to
  # fields in the original cell. For example, new.lattice
  # points to old.lattice. A deepcopy needs to be done to 
  # untangle the two, but I'm unsure if I care to do it.
  Cell(
    cell.lattice, spos, cell.velocity, cell.masses, 
    cell.symbols, cell.mask, cell.PBC, cell.NC
  )
end

"""
    HiddenOptVars

Structure for hidden variable optimization.

# Fields
- `potVars`: Potential variables.
- `cellBuf`: Cell buffer.
- `superBuf`: Supercell buffer.
- `scaleEnergy`: Energy scaling factor.
- `PBC`: Periodic boundary conditions.
- `mols`: Molecule indices.
- `pars`: Pair indices.
- `T`: Transformation matrix.
"""
struct HiddenOptVars
  potVars::PotVars
  cellBuf::MyCell
  superBuf::MyCell
  scaleEnergy::Float64
  PBC::Vector{Bool}
  mols::Vector
  pars::Vector
  T::Matrix
end

"""
    hiddenOpt(calc::MyCalc, algo, cell, T; kwargs...)

Optimize a cell using a hidden variable approach and a supercell transformation.

# Arguments
- `calc`: Calculator object.
- `algo`: Optimization algorithm.
- `cell`: `MyCell` object.
- `T`: Transformation matrix.
- `kwargs`: Additional options.

# Returns
- Optimized `MyCell` object.
"""
function hiddenOpt(calc::MyCalc, algo, cell, T; kwargs...)

  super   = makeSuperCell(cell, T)
  y, vars = prep4pot(calc.b, super)
  x0      = prepX0(cell)
  Γ       = zero(y)

  hideVars = HiddenOptVars(
    vars.potVars,
    deepcopy(cell),
    super,
    length(cell.masses) / length(super.masses),
    cell.PBC,
    vars.mols,
    vars.pars,
    T
  )

  optFunc  = NLSolversBase.only_fg!((F,G,x) -> hiddenEoM(F,G,Γ,calc,hideVars,x))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  spos     = getScaledPos(res.minimizer, cell.lattice)
  spos     = [MVector(i) for i in spos]
  

  # The new cell returned will have fields that point to
  # fields in the original cell. For example, new.lattice
  # points to old.lattice. A deepcopy needs to be done to 
  # untangle the two, but I'm unsure if I care to do it.
  Cell(
    cell.lattice, spos, cell.velocity, cell.masses,
    cell.symbols, cell.mask, cell.PBC, cell.NC
  )  
end

"""
    hiddenEoM(F, G, Γ, calc::MyCalc, vars, x)

Energy and gradient evaluation for hidden variable optimization.

# Arguments
- `F`: Energy output (or nothing).
- `G`: Gradient output (or nothing).
- `Γ`: Gradient buffer.
- `calc`: Calculator object.
- `vars`: HiddenOptVars object.
- `x`: Coordinate vector.

# Returns
- Energy value if `F` is not `nothing`.
"""
function hiddenEoM(F, G, Γ, calc::MyCalc, vars, x)
  getScaledPos!(vars.cellBuf, x)

  makeSuperCell!(vars.superBuf, vars.cellBuf, vars.T)

  y = prepX0(vars.superBuf)
  E = fg!(F, Γ, y, vars, calc)

  if G != nothing
    G .= Γ[1:length(G)]
  end
  
  if F != nothing
    return E * vars.scaleEnergy
  end

end

"""
    optCell(calc::MyCalc, algo, cell::MyCell; precon=nothing, kwargs...)

Optimize the lattice parameters of a cell.

# Arguments
- `calc`: Calculator object.
- `algo`: Optimization algorithm.
- `cell`: `MyCell` object.
- `precon`: (Optional) Preconditioner.
- `kwargs`: Additional options.

# Returns
- Optimized `MyCell` object.
"""
function optCell(calc::MyCalc, algo, cell::MyCell; precon=nothing, kwargs...)

  ret          = deepcopy(cell)
  diag         = isdiag(ret.lattice)
  lat0         = reshape(ret.lattice, 9) |> Vector
  optFunc      = NLSolversBase.only_fg!(
    (F,G,x) -> st!(F,G,x, ret, calc; precon=precon, diag=diag)
  )
  convCrit     = Optim.Options(; kwargs...)
  res          = optimize(optFunc, lat0, algo, convCrit)
  newLat       = reshape(res.minimizer, (3,3))
  ret.lattice .= newLat

  ret
end