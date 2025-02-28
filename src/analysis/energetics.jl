function vibExcite!(mol, eignvec, E)
  M = [i.m for i in mol for j in 1:3]
  v = @. sqrt( 2E / M ) * eignvec

  vcom = vCoM(mol)
  for i in 1:3:length(v)
    j::UInt32 = (i+2)/3
    mol[j].v += v[i:i+2]
  end 

  if !isapprox(vCoM(mol), vcom; atol=1e-8)
    println("Uh Oh")
    zeroVCoM!(mol)
  end

end

function transExcite!(mol, KE)
  r  = randVector()
  ke = KE * 8.6173e-5 #convert Kelvin to eV

  for i in mol
    i.v += sqrt(2ke / i.m) .* r
  end
end

function getTransEnergy(mol)
  μ = reducedMass(mol)
  v = vCoM(mol)
  E = 0.5 * μ * dot(v,v)
  
  E
end

function getTransEnergy(pos,vel,mas)
  μ = reducedMass(mas)
  v = vCoM(vel,mas)
  E = 0.5 * μ * dot(v,v)
  
  E
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

  E
end

function getRotEnergy(pos,vel,mas)
  vcom = vCoM(vel,mas)
  com  =  CoM(pos,mas)
  E    = 0.0

  for i in 1:length(mas)
    r  = pos[i] - com
    v  = vel[i] - vcom
    w  = cross(r,v) / dot(r,r)
    I  = mas[i] * dot(r,r)
    E += 0.5 * I * dot(w,w)
  end
  
  E
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

  E
end

#The general method above will be great for more complicated molecules
# but it is very slow. While below only works for CO but is fast.

function getCOVibEnergy(mol; pot=nothing)
  diff = mol[2].r - mol[1].r
  rhat = normalize(diff)
  v1   = dot(mol[1].v, rhat)
  v2   = dot(mol[2].v, rhat)
  E    = 0.5*mol[1].m*v1^2 + 0.5*mol[2].m*v2^2

  if pot != nothing
    E += getPotEnergy(pot, mol)
  end

  E
end

function getCOVibEnergy(pos,vel,mas; pot=nothing)
  diff = pos[2] - pos[1]
  rhat = normalize(diff)
  v1   = dot(vel[1], rhat)
  v2   = dot(vel[2], rhat)
  E    = 0.5*mas[1]*v1^2 + 0.5*mas[2]*v2^2

  if pot != nothing
    x0 = [j for i in pos for j in i]
    E += pot(true, nothing, x0, optVars([[1,2]], []))
  end

  E
end

function getPotEnergy(EoM, bdys::Vector{MyAtoms})
  x0, vars = prep4pot(EoM, bdys)
  energy   = EoM(true, nothing, x0, vars)
  
  energy
end

function getPotEnergy(EoM, cell::MyCell)
  x0, vars = prep4pot(EoM, cell)
  energy   = EoM(true, nothing, x0, vars)
  
  energy
end