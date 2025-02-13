
struct optVars
  potVars::PotVars
  mols::Vector
  pars::Vector
  m::Vector
  PBC::Vector{Bool}
  NC::Vector{Int32}
  lattice::Matrix
end


function prepX0(bdys)
  r  = [i.r for i in bdys]
  x0 = [j for i in r for j in i]
  
  x0
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

#TODO: Update this to stop making bdys
function prep4pot(EoM, cell::MyCell)
  bdys       = makeBdys(cell)
  x0         = prepX0(bdys)
  potVars    = EoM(cell)
  pars, mols = getPairs(bdys)
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

  Cell(cell.lattice, spos, cell.masses, cell.symbols, cell.PBC, cell.NC)
end

function optCell(EoM, cell; kwargs...)

  # Need to init vars for simulations 
  p    = makeBdys(cell) |> EoM
  #TODO: make EoM(cell) variant to avoid making bodies above

  x0       = cell.scaled_pos
  # box0     = transpose(cell.lattice) |> (x -> reshape(x, 9)) |> (
  #   x -> convert(Vector{Float64}, x))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(c -> EoM(cell.lattice .* c, x0, p), [1.0], NelderMead(), convCrit)
  newBox   = cell.lattice .* res.minimizer

  Cell(newBox, cell.scaled_pos, cell.masses, cell.symbols, cell.PBC, cell.NC)
end