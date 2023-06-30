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
  # Method taken from https://mathworld.wolfram.com/SpherePointPicking.html
  # only I swap phi and theta 

  u = rand()
  v = rand()
  θ = acos(2v - 1)
  ϕ = 2pi*u
  
  x = cos(ϕ) * sin(θ)
  y = sin(ϕ) * sin(θ)
  z = cos(θ)

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

function spawnMol(mol, bdys, com)
  v     = randVector()
  d     = maximum([dot(i.r-com, v) for i in bdys]) + 8
  R     = d * v + com
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

function hitAndStick(EoM, inp; callback=nothing)

  #Initialize System at Center
  bdys = readXyz(inp.xyz)
  mol  = readXyz(inp.mol)

  while length(bdys) < inp.size
    #Get current center of mass
    com = CoM(bdys)

    #Randomly rotate incoming molecule
    randRotate!(mol) 
    #If i dont reset mol, then the rotations are cummulative
    # is this more random??

    #Spawn New Molecule
    new = spawnMol(mol, bdys, com)

    #Update velocities
    giveKE!(new, com, inp.KE)

    #Add new to bdys
    push!(bdys, new...)
    
    #Run NVE
    time = inp.htime * ps
    solu = runNVE(EoM, (0, time), fs, bdys)
    getLastFrame!(bdys, solu)
    
    #Free memory
    solu = 0
    GC.gc()

    #Run NVT
    time = inp.stime * ps
    solu = runNVT(EoM, (0, time), fs, bdys, inp.thermo, inp.thermoInps)
    getLastFrame!(bdys, solu)
    
    #Free memory
    solu = 0
    GC.gc()

    #If needed, execute callback function
    if callback != nothing
      callback(bdys)
    end

  end

  writeXyz(inp.save, bdys)
end
