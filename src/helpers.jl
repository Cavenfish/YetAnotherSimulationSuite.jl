
function getLastFrame(solu)
  n = length(solu.prob.p.bdys)
  
  new = Atom[]
  for i in 1:n
    r = solu.u[end].x[2][i]
    v = solu.u[end].x[1][i]
    m = solu.prob.p.bdys[i].m
    s = solu.prob.p.bdys[i].s
    push!(new, Atom(r, v, m, s))
  end
  return new
end

function getLastFrame!(bdys, solu)
  N = length(bdys)
  
  for i in 1:N
    bdys[i].r = solu.u[end].x[2][i]
    bdys[i].v = solu.u[end].x[1][i]
    bdys[i].m = solu.prob.p.bdys[i].m
    bdys[i].s = solu.prob.p.bdys[i].s
  end
end

function CoM(bdys)
  M = sum([i.m for i in bdys])
  r = sum([i.m*i.r for i in bdys])
  return r ./ M
end

function reducedMass(bdys)
  u = sum([i.m^-1 for i in bdys])^-1
  return u
end

function vibExcite!(mol, eignvec, E)
  v = sqrt.((2*E) ./ mol.m) # column vector of length of atoms in mol

  mol.v .+= v .* eignvec
end
