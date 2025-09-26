struct HnS
  size::UInt16
  htime::Quantity
  stime::Quantity
  KE::Quantity
  xyz::String
  mol::String
  save::String
  dt::Quantity
  thermo::MyThermostat
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
  
  r
end

function spawnMol(mol, bdys, com)
  v     = randVector()
  d     = maximum([dot(i.r-com, v) for i in bdys]) + 8
  R     = d * v + com
  spawn = deepcopy(mol)

  translate!(spawn, R)

  spawn
end

function giveKE!(mol, com, KE)
  r  = com - CoM(mol)
  r /= norm(r)
  ke = uconvert(u"eV", KE) |> ustrip

  for i in mol
    i.v += sqrt(2ke / i.m) .* r
  end
end

function hitAndStick(EoM, inp; callback=nothing)

  nve = NVE()
  nvt = inp.thermo |> NVT

  #Initialize System at Center
  bdys = readSystem(inp.xyz)
  mol  = readSystem(inp.mol)

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
    solu = run(EoM, bdys, inp.htime, inp.dt, nve)
    getLastFrame!(bdys, solu)
    zeroVCoM!(bdys)
    
    #Free memory
    @free solu

    #Run NVT
    solu = run(EoM, bdys, inp.htime, inp.dt, nvt)
    getLastFrame!(bdys, solu)
    
    #Free memory
    @free solu

    #If needed, execute callback function
    if callback != nothing
      callback(bdys)
    end

  end

  write(inp.save, bdys)
end
