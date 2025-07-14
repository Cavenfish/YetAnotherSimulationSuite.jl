# TODO:
#   - Generalize the dimensionality of optimizations

struct optVars{D,B,P, I<:Int, PV<:PotVars, F<:AbstractFloat}
  potVars::PV
  mols::Vector{Vector{I}}
  pars::Vector{P}
  m::Vector{F}
  PBC::Vector{B}
  NC::Vector{I}
  lattice::MMatrix{D,D,F}
end


function prepX0(bdys::Vector{MyAtoms})
  r  = [i.r for i in bdys]
  
  [j for i in r for j in i]
end

function prepX0(cell::MyCell)
  r  = getPos(cell)
  
  [j for i in r for j in i]
end

function prep4pot(builder, bdys::Vector{MyAtoms})
  m          = [i.m for i in bdys]
  x0         = prepX0(bdys)
  potVars    = builder(bdys)
  pars, mols = getPairs(bdys)
  NC         = [0,0,0]
  PBC        = repeat([false], 3)
  lattice    = MMatrix{3,3}(zeros(3,3))
  vars       = optVars(potVars, mols, pars, m, PBC, NC, lattice)
  
  x0, vars
end

function prep4pot(builder, cell::MyCell)
  bdys       = makeBdys(cell)
  x0         = prepX0(cell)
  potVars    = builder(cell)
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

    r = MVector{3}(opt[i:i+2])
    v = bdys[j].v
    m = bdys[j].m
    s = bdys[j].s
    push!(new, Particle(r,v,m,s))
  end

  new
end

function opt(calc::MyCalc, algo, bdys::Vector{MyAtoms}; kwargs...)
  x0, vars = prep4pot(calc.b, bdys)
  optFunc  = Optim.only_fg!((F,G,x) -> fg!(F,G,x, vars, calc))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  optBdys  = getNewBdys(bdys, res)

  optBdys
end

function opt(calc::MyCalc, algo, cell::MyCell; kwargs...)
  x0, vars = prep4pot(calc.b, cell)
  optFunc  = Optim.only_fg!((F,G,x) -> fg!(F,G,x, vars, calc))
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

  optFunc  = Optim.only_fg!((F,G,x) -> hiddenEoM(F,G,Γ,calc,hideVars,x))
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

function optCell(calc::MyCalc, algo, cell::MyCell; precon=nothing, kwargs...)

  ret          = deepcopy(cell)
  diag         = isdiag(cell.lattice)
  lat0         = reshape(cell.lattice, 9) |> Vector
  optFunc      = Optim.only_fg!(
    (F,G,x) -> st!(F,G,x, cell, calc; precon=precon, diag=diag)
  )
  convCrit     = Optim.Options(; kwargs...)
  res          = optimize(optFunc, lat0, algo, convCrit)
  newLat       = reshape(res.minimizer, (3,3))
  ret.lattice .= newLat

  ret
end