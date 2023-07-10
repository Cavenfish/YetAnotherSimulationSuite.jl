
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

function vCoM(bdys)
  M = sum([i.m for i in bdys])
  v = sum([i.m*i.v for i in bdys])
  return v ./ M
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

function getTransEnergy(mol)
  μ = reducedMass(mol)
  v = vCoM(mol)
  E = 0.5 * μ * dot(v,v)
  return E
end

function getRotEnergy(mol)
  vcom = vCoM(mol)
  com  =  CoM(mol)
  E    = 0.0

  for i in mol
    r  = i.r - com
    v  = i.v - vcom
    w  = cross(r,v) / dot(r,r)
    I  = i.m * dot(r,r) 
    E += 0.5 * I * dot(w,w)
  end

  return E
end

function getVibEnergy(mol, eignvec; pot=nothing)

  E = 0.0
  for i in 1:3:length(eignvec)
    j::UInt32 = (i+2)/3
    ehat      = eignvec[i:i+2] / norm(eignvec[i:i+2])
    E        += 0.5 * mol[j].m * dot(mol[j].v, ehat)^2
  end

  if pot != nothing
    E += getPotEnergy(pot, mol)
  end

  return E
end

#The general method above will be great for more complicated molecules
# but it is very slow. While below only works for CO but is fast.

function getCOVibEnergy(mol; pot=nothing)
  diff = mol[2].r - mol[1].r
  rhat = diff / norm(diff)
  v1   = dot(mol[1].v, rhat)
  v2   = dot(mol[2].v, rhat)
  E    = 0.5*mol[1].m*v1^2 + 0.5*mol[2].m*v2^2

  if pot != nothing
    E += getPotEnergy(pot, mol)
  end

  return E
end

function getPotEnergy(EoM, bdys)
  x0, vars = prep4pot(bdys)
  energy   = EoM(true, nothing, x0, vars)
  return energy
end

