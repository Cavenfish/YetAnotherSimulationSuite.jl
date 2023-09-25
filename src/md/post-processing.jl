
struct Traj <: MyTraj
  t::Vector{Float64}
  E::Vector{Float64}
  T::Vector{Float64}
  F::Vector
  r::Vector
  v::Vector
  m::Vector{Float64}
end

function processDynamics(solu; dt=fs, step=1)
  t = solu.t[1:step:end] / dt
  N = length(solu.t)

  if typeof(solu.prob.p.energy[1]) == NTuple{5, Float64}
    E = [i[1] for i in solu.prob.p.energy[1:step:end]]
  else
    E = solu.prob.p.energy[1:step:end]
  end

  if typeof(solu.prob.p.forces[1]) == NTuple{5, Float64}
    F = [i[1] for i in solu.prob.p.forces[1:step:end]]
  else
    F = solu.prob.p.forces[1:step:end]
  end

  r = [solu.u[i].x[2] for i in 1:step:N]
  v = [solu.u[i].x[1] for i in 1:step:N]
  m = solu.prob.p.m
  T = [getTemp(m, i, kB, length(m)) for i in v]

  return Traj(t, E, T, F, r, v, m)
end

function processDynamics!(tj, solu; dt=fs, step=1)
  t = solu.t[1:step:end] / dt
  N = length(solu.t)

  if typeof(solu.prob.p.energy[1]) == NTuple{5, Float64}
    E = [i[1] for i in solu.prob.p.energy[1:step:end]]
  else
    E = solu.prob.p.energy[1:step:end]
  end

  if typeof(solu.prob.p.forces[1]) == NTuple{5, Float64}
    F = [i[1] for i in solu.prob.p.forces[1:step:end]]
  else
    F = solu.prob.p.forces[1:step:end]
  end

  r = [solu.u[i].x[2] for i in 1:step:N]
  v = [solu.u[i].x[1] for i in 1:step:N]
  m = solu.prob.p.m
  T = [getTemp(m, i, kB, length(m)) for i in v]

  push!(tj.t, t...)
  push!(tj.E, E...)
  push!(tj.T, T...)
  push!(tj.F, F...)
  push!(tj.r, r...)
  push!(tj.v, v...)
end

function processTmpFiles(files; kwargs...)
  
  f    = open(files[1], "r")
  solu = deserialize(f)
  tj   = processDynamics(solu; kwargs...) 
  close(f)

  # Free memory
  @free solu

  for file in files[2:end]
    f    = open(file, "r")
    solu = deserialize(f)
    processDynamics!(tj, solu; kwargs...)
    close(f)

    #Free memory
    @free solu
  end

  return tj
end

function trackVACF(files, safe)

  df = DataFrame()

  for file in files
    f    = open(file, "r")
    solu = deserialize(f)
    close(f)

    v,m  = getVelMas(solu)

    inp  = vacfInps(v[1:safe], m, 1e15, false, HannM, 8, true)
    out  = VDOS(inp)
    col  = replace(file, ".tmp" => "")

    if isempty(df)
      df[!, "v"] = out.v
    end

    df[!, col] = out.I

    @free solu
  end

  return df
end

function trackEnergyDissipation(traj, pot, mol)

  m = traj.m
  N = length(traj.t)
  n = length(m)

  molVib, molRot, molTra = Float64[], Float64[], Float64[]
  avgVib, avgRot, avgTra = Float64[], Float64[], Float64[]

  others = [[i,i+1] for i in 1:2:n if !(i in mol)]

  for i in 1:N
    r = traj.r[i]
    v = traj.v[i]

    mvib = getCOVibEnergy(r[mol], v[mol], m[mol]; pot=pot)
    mrot = getRotEnergy(  r[mol], v[mol], m[mol])
    mtra = getTransEnergy(r[mol], v[mol], m[mol])

    push!(molVib, mvib)
    push!(molRot, mrot)
    push!(molTra, mtra)

    avib,arot,atra = 0.0,0.0,0.0
    for j in others
      avib += getCOVibEnergy(r[j], v[j], m[j]; pot=pot)
      arot += getRotEnergy(  r[j], v[j], m[j])
      atra += getTransEnergy(r[j], v[j], m[j])
    end


    push!(avgVib, avib/n)
    push!(avgRot, arot/n)
    push!(avgTra, atra/n)
  end

  df = DataFrame(  time=traj.t,
                   temp=traj.T,
                 molVib=molVib,
                 molRot=molRot,
                 molTra=molTra,
                 avgVib=avgVib,
                 avgRot=avgRot,
                 avgTra=avgTra)

  return df
end

function trackAllVibEnergy(tj; pot=nothing)

  N   = length(tj.r[1])

  pos = tj.r[102]
  vel = tj.v[102]
  mas = tj.m

  tmp = [getCOVibEnergy(pos[i:i+1], vel[i:i+1], mas[i:i+1]) for i in 1:2:N]

  # This is the excited molecule
  ind = findall(e -> e == maximum(tmp), tmp)[1]

  dis = DataFrame()
  vib = DataFrame()

  dis[!, "time"] = tj.t
  vib[!, "time"] = tj.t

  for i in 1:div(N,2)
    str = "mol$i"
    
    dis[!, str] = zero(tj.t)
    vib[!, str] = zero(tj.t)
  end

  for i in 1:length(tj.t)
    for j in 1:div(N,2)
      #properties of excited molecule
      epos = tj.r[i][ind*2-1:ind*2]
      emas = tj.m[ind*2-1:ind*2]
      ecom = CoM(epos, emas)

      a   = j*2-1
      b   = j*2

      #properties of other molecule
      pos = tj.r[i][a:b]
      vel = tj.v[i][a:b]
      mas = tj.m[a:b]
      com = CoM(pos, mas)
      
      #distance and vib energy
      d = norm(com - ecom)
      v = getCOVibEnergy(pos, vel, mas; pot=pot)

      #store in DFs
      dis[i,j+1] = d
      vib[i,j+1] = v

    end
  end

  return dis, vib
end

function trackRadialEnergy(tj; pot=nothing)

  #TODO: Make this function collect vib, tra, and rot energy at radial slices

  
  
end