
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

  for i in 1:N
    a = 4*i - 3
    b = 3*i
    r = SVector{3}(opt[a:b])
    v = bdys[i].v
    m = bdys[i].m
    s = bdys[i].s
    push!(new, Atom(r,v,m,s))
  end

  return new
end

function opt(EoM, algo, bdys)

  x0         = prep_x0(bdys)
  pars, mols = getPairs(bdys)
  vars       = optVars(mols, pars)

  res        = optimize(x -> EoM(x, vars), x0, algo)
  optBdys    = getNewBdys(bdys, res)

  return optBdys
end
