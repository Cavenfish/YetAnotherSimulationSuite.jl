struct HnS
  size::UInt16
  htime::UInt16
  stime::UInt16
  KE::Float64
  xyz::String
  mol::String
  save::String
  thermo
  thermoInps
end

function randVector()
  R = 1
  θ = rand() * pi
  ϕ = rand() * 2 * pi 
  
  x = R * cos(ϕ) * sin(θ)
  y = R * sin(ϕ) * sin(θ)
  z = R * cos(θ)

  r  = [x,y,z]
  r /= norm(r) 
  return r
end

function randRotate!(mol)
  α = rand(-pi:1e-10:pi)
  γ = rand(-pi:1e-10:pi)
  β  = rand(0:1e-10:pi)

  Rz = [cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
  Ry = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)]
  Rx = [1 0 0; 0 cos(γ) -sin(γ); 0 sin(γ) cos(γ)]

  R  = Rz*Ry*Rx

  for i in mol
    i.r = R * i.r
  end
end

function spawnMol(mol, d, com)
  R     = d * randVector() + com
  spawn = Atom[]

  for i in mol
    r = i.r + R
    push!(spawn, Atom(r, i.v, i.m, i.s))
  end

  return spawn
end

function giveKE!(mol, com, KE)
  r  = com - CoM(mol)
  r /= norm(r)
  ke = KE * 8.6173e-5 #convert Kelvin to eV

  for i in mol
    i.v += sqrt(2ke / i.m) .* r
  end
end

function hitAndStick(EoM, inp)

  #Initialize System at Center
  bdys = readXyz(inp.xyz)
  mol  = readXyz(inp.mol)

  while length(bdys) < inp.size
    #Get current center of mass
    com = CoM(bdys)

    #Get furthest molecule from CoM
    d = maximum([norm(i.r-com) for i in bdys]) + 8

    #Randomly rotate incoming molecule
    randRotate!(mol) 
    #If i dont reset mol, then the rotations are cummulative
    # is this more random??

    #Spawn New Molecule
    new = spawnMol(mol, d, com)

    #Update velocities
    giveKE!(new, com, inp.KE)

    #Add new to bdys
    push!(bdys, new...)
    
    #Run NVE
    time = inp.htime * ps
    solu = runNVE(EoM, (0, time), fs, bdys)
    bdys = getLastFrame(solu)
    
    #Run NVT
    time = inp.stime * ps
    solu = runNVT(EoM, (0, time), fs, bdys, inp.thermo, inp.thermoInps)
    bdys = getLastFrame(solu)

  end

  writeXyz(inp.save, bdys)
end
