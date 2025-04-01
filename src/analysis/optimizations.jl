
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

  Cell(cell.lattice, spos, cell.velocity, cell.masses, cell.symbols, cell.PBC, cell.NC)
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


# optSuper and potPBC are good first run functions 
# however, they are not efficient and NEED to be redone.

# Should be done without reallocating cells/potVars at
# every iteration.

function optSuper(EoM, algo, cell, T; kwargs...)

  x0       = getPos(cell) |> (x -> vcat(x...))
  tmp      = deepcopy(cell)
  optFunc  = Optim.only_fg!((F,G,x) -> potPBC(F,G, EoM, tmp, T, x))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  spos     = getScaledPos(res.minimizer, cell.lattice)

  Cell(cell.lattice, spos, cell.velocity, cell.masses, cell.symbols, cell.PBC, cell.NC)  
end

function potPBC(F, G, EoM, cell, T, x0)
  spos = getScaledPos(x0, cell.lattice)

  for i in 1:length(cell.scaled_pos)
    cell.scaled_pos[i] .= spos[i]
  end

  N     = length(x0)
  super = makeSuperCell(cell, T)

  if G != nothing
    G .= - getForces(EoM, super) |> (x -> vcat(x...)[1:N])
  end
  
  if F != nothing
    return getPotEnergy(EoM, super) / det(T)
  end

end