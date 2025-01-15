
struct optVars
  potVars::PotVars
  mols::Vector
  pars::Vector
  m::Vector
end


function prepX0(bdys)
  r  = [i.r for i in bdys]
  x0 = [j for i in r for j in i]
  
  x0
end

function prep4pot(EoM, bdys)
  m          = [i.m for i in bdys]
  x0         = prepX0(bdys)
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  vars       = optVars(potVars, mols, pars, m)
  
  x0, vars
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

  new
end

function getNewBdys(bdys, x0, box)
  N    = length(bdys)
  new  = Atom[]
  x0 .*= repeat(box, N)

  for i in 1:3:N*3
    j::UInt16 = (i+2) / 3

    r = SVector{3}([x0[i], x0[i+1], x0[i+2]])
    v = bdys[j].v
    m = bdys[j].m
    s = bdys[j].s
    push!(new, Atom(r,v,m,s))
  end

  new
end

function opt(EoM, algo, bdys; kwargs...)

  x0, vars = prep4pot(EoM, bdys)
  optFunc  = Optim.only_fg!((F,G,x) -> EoM(F,G,x, vars))
  convCrit = Optim.Options(; kwargs...)
  res      = optimize(optFunc, x0, algo, convCrit)
  optBdys  = getNewBdys(bdys, res)

  optBdys
end

function optCell(EoM, bdys, box0; kwargs...)

  x0, vars = prep4pot(EoM, bdys)

  for i = 1:3:length(x0)
    x0[i]   /= box0[1]
    x0[i+1] /= box0[2]
    x0[i+2] /= box0[3]
  end

  convCrit = Optim.Options(; kwargs...)
  res      = optimize(box -> EoM(box, x0), box0, NelderMead(), convCrit)
  newBox   = res.minimizer
  newBdys  = getNewBdys(bdys, x0, newBox)

  newBdys, newBox
end