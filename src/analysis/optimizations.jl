
struct optVars
  potVars::PotVars
  mols::Vector
  pars::Vector
  m::Vector
  PBC::Vector{Bool}
  NC::Vector{Int32}
  lattice::Matrix
end


function prepX0(bdys::Vector{MyAtoms})
  r  = [i.r for i in bdys]
  
  [j for i in r for j in i]
end

function prepX0(cell::MyCell)
  r  = getPos(cell)
  
  [j for i in r for j in i]
end

function prep4pot(EoM, bdys::Vector{MyAtoms})
  m          = [i.m for i in bdys]
  x0         = prepX0(bdys)
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  NC         = [0,0,0]
  PBC        = repeat([false], 3)
  lattice    = zeros(3,3)
  vars       = optVars(potVars, mols, pars, m, PBC, NC, lattice)
  
  x0, vars
end

function prep4pot(EoM, cell::MyCell)
  bdys       = makeBdys(cell)
  x0         = prepX0(cell)
  potVars    = EoM(cell)
  pars, mols = getPairs(cell)
  vars       = optVars(potVars, mols, pars, cell.masses, 
                       cell.PBC, cell.NC, cell.lattice)
  
  x0, vars
end

function getNewBdys(bdys, res)
  N   = length(bdys)
  opt = res.minimizer
  new = MyAtoms[]

  for i in 1:3:N*3
    j::UInt16 = (i+2) / 3

    r = SVector{3}(opt[i:i+2])
    v = bdys[j].v
    m = bdys[j].m
    s = bdys[j].s
    push!(new, Atom(r,v,m,s))
  end

  new
end

function opt(EoM, algo, bdys::Vector{MyAtoms}; kwargs...)

  x0, vars = prep4pot(EoM, bdys)
  optFunc  = Optim.only_fg!((F,G,x) -> EoM(F,G,x, vars))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  optBdys  = getNewBdys(bdys, res)

  optBdys
end

function opt(EoM, algo, cell::MyCell; kwargs...)

  x0, vars = prep4pot(EoM, cell)
  optFunc  = Optim.only_fg!((F,G,x) -> EoM(F,G,x, vars))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  spos     = getScaledPos(res.minimizer, cell.lattice)

  # The new cell returned will have fields that point to
  # fields in the original cell. For example, new.lattice
  # points to old.lattice. A deepcopy needs to be done to 
  # untangle the two, but I'm unsure if I care to do it.
  Cell(
    cell.lattice, spos, cell.velocity, cell.masses, 
    cell.symbols, cell.mask, cell.PBC, cell.NC
  )
end

function optCell(EoM, algo, cell::MyCell; kwargs...)

  lat0         = reshape(cell.lattice, 9)
  optFunc      = Optim.only_fg!((F,G,x) -> EoM(F,G, cell, x))
  convCrit     = Optim.Options(; kwargs...)
  res          = optimize(optFunc, lat0, algo, convCrit)
  newLat       = reshape(res.minimizer, (3,3))
  ret          = deepcopy(cell)
  ret.lattice .= newLat

  ret
end

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

function hiddenOpt(EoM, algo, cell, T; kwargs...)

  super   = makeSuperCell(cell, T)
  y, vars = prep4pot(EoM, super)
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

  optFunc  = Optim.only_fg!((F,G,x) -> hiddenEoM(F,G,Γ,EoM,hideVars,x))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  spos     = getScaledPos(res.minimizer, cell.lattice)

  # The new cell returned will have fields that point to
  # fields in the original cell. For example, new.lattice
  # points to old.lattice. A deepcopy needs to be done to 
  # untangle the two, but I'm unsure if I care to do it.
  Cell(
    cell.lattice, spos, cell.velocity, cell.masses,
    cell.symbols, cell.mask, cell.PBC, cell.NC
  )  
end

function hiddenEoM(F, G, Γ, EoM, vars, x)
  getScaledPos!(vars.cellBuf, x)

  makeSuperCell!(vars.superBuf, vars.cellBuf, vars.T)

  y = prepX0(vars.superBuf)
  E = EoM(F, Γ, y, vars)

  if G != nothing
    G .= Γ[1:length(G)]
  end
  
  if F != nothing
    return E * vars.scaleEnergy
  end

end