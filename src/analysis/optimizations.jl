
struct optVars
  mols::Vector
  pars::Vector
end


function prep_x0(bdys)
  r  = [i.r for i in bdys]
  x0 = [j for i in r for j in i]
  return x0
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

  x0         = prep_x0(bdys)
  pars, mols = getPairs(bdys)
  vars       = optVars(mols, pars)
  optFunc    = Optim.only_fg!((F,G,x) -> EoM(F,G,x, vars))
  res        = optimize(optFunc, x0, algo; kwargs...)
  optBdys    = getNewBdys(bdys, res)

  return optBdys
end