
struct optVars
  mols::Vector
  pars::Vector
  m::Vector
end


function prepX0(bdys)
  r  = [i.r for i in bdys]
  x0 = [j for i in r for j in i]
  return x0
end

function prep4pot(bdys)
  m          = [i.m for i in bdys]
  x0         = prepX0(bdys)
  pars, mols = getPairs(bdys)
  vars       = optVars(mols, pars, m)
  return x0, vars
end

function getNewBdys(bdys, res)
  N   = length(bdys)
  opt = res.minimizer
  new = Atom[]

  for i in 1:3:N*3
    j::UInt16 = (i+2) / 3

    r = SVector{3}(opt[i:i+2])
    v = bdys[j].v
    m = bdys[j].m
    s = bdys[j].s
    push!(new, Atom(r,v,m,s))
  end

  return new
end

function opt(EoM, algo, bdys; kwargs...)

  m          = [i.m for i in bdys]
  x0         = prepX0(bdys)
  pars, mols = getPairs(bdys)
  vars       = optVars(mols, pars, m)
  optFunc    = Optim.only_fg!((F,G,x) -> EoM(F,G,x, vars))
  convCrit   = Optim.Options(; kwargs...)
  res        = optimize(optFunc, x0, algo, convCrit)
  optBdys    = getNewBdys(bdys, res)

  return optBdys
end