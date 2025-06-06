function getPotEnergy(calc::MyCalc, obj::Union{MyCell, Vector{MyAtoms}})
  x, vars = prep4pot(calc.b, obj)
  energy  = fg!(true, nothing, x, vars, calc)
  
  energy
end

function getVibEnergy(mol::Vector{MyAtoms}, eignvec; calc=nothing)

  E = 0.0
  for i in 1:3:length(eignvec)
    j::UInt32 = (i+2)/3
    ehat      = eignvec[i:i+2] / norm(eignvec[i:i+2])
    E        += 0.5 * mol[j].m * dot(mol[j].v, ehat)^2
  end

  if calc != nothing
    E += getPotEnergy(calc, mol)
  end

  E
end

function getTransEnergy(mol::Vector{MyAtoms})
  μ = reducedMass(mol)
  v = vCoM(mol)
  E = 0.5 * μ * dot(v,v)
  
  E
end

function getRotEnergy(mol::Vector{MyAtoms})
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

function vibExcite!(mol::Vector{MyAtoms}, eignvec, E)
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

function transExcite!(mol::Vector{MyAtoms}, ke)
  r = randVector()
  μ = reducedMass(mol)
  v = sqrt(2ke / μ)

  for i in mol
    i.v .+= v .* r
  end
end