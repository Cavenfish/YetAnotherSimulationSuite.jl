
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

function swapIso!(bdys, swap, mas)
  for i in swap
    bdys[i].m = mas[i]
  end
end

function vibExcite!(mol, eignvec, E)
  M = [i.m for i in mol for j in 1:3]
  v = @. sqrt( 2E / M ) * eignvec

  for i in 1:3:length(v)
    j::UInt32 = (i+2)/3
    mol[j].v += v[i:i+2]
  end 

end

function getVibEnergy(mol, eignvec)
  v = [i.v[j] for i in mol for j in 1:3]
  m = [i.m    for i in mol for j in 1:3]

  E = 0.5 * dot((m .^ 0.5) .* v, eignvec)^2
  return E
end

function getPotEnergy(EoM, bdys)
  x0, vars = prep4pot(bdys)
  energy   = EoM(true, nothing, x0, vars)
  return energy
end

