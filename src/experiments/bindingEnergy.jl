"""
Calculates binding energies based on given
toml input card.

Example of Input Card
-----------------------

[Settings]
EoM = "COCO"
mol = "/home/brian/COjl/xyzFiles/1.xyz"
clu = "/home/brian/COjl/clu.xyz"
sites = 250
saveat = 10
minD = 3.5

[OPT]
algo = "LBFGS"
pre.x_tol = 1e-5
pre.f_tol = 1e-12
pre.g_tol = 1e-6
pre.iterations = 100000
pst.x_tol = 1e-5
pst.f_tol = 1e-12
pst.g_tol = 1e-3
pst.iterations = 100000

[Saving]
df = "beDF.jld2"
"""

function calcBEs(inpFile::String)
  inp = TOML.parsefile(inpFile)

  #Split dict for easier usage
  cnfg = inp["Settings"]

  #Load up EoM and opt algo
  EoM  = mkvar(cnfg["EoM"])
  algo = mkvar(inp["OPT"]["algo"])()

  #Build pre opt kwargs dict
  pre = Dict()
  for key in keys(inp["OPT"]["pre"])
    value = inp["OPT"]["pre"][key]
    pre[Symbol(key)] = value
  end

  #Build post opt kwargs dict
  pst = Dict()
  for key in keys(inp["OPT"]["pst"])
    value = inp["OPT"]["pst"][key]
    pst[Symbol(key)] = value
  end

  #Get leftover vars
  minD  = cnfg["minD"]
  N     = cnfg["sites"]
  savN  = cnfg["saveat"]
  nthrd = Threads.nthreads()

  #Read in and pre-optimise cluster and mol
  mol = readXyz(cnfg["mol"]) |> (x -> opt(EoM, algo, x; pre...))
  clu = readXyz(cnfg["clu"])  |> (x -> opt(EoM, algo, x; pre...))

  println("pre-opt done")

  #Get cluster center of mass (for surface norms later)
  com = CoM(clu)

  # Get energies
  molE = getPotEnergy(EoM, mol)
  cluE = getPotEnergy(EoM, clu)

  #Get spots on surface of cluster
  spots = [i.r for i in clu] |> alphashape |> getSpots |> (x -> sample(x, N; replace=false))

  #Get rough count of iterations per thread with padding for faster threads
  tmp = div(savN, nthrd) + 2

  #Check for threading efficiency
  savN <= nthrd || @warn "SaveAt is less than nthreads.\n This is not efficient!"

  #Prep BE array
  ret = zeros(nthrd, tmp)
  BE  = [0.0 for i in 1:N]
  i   = 1

  while i < N

    Threads.@threads for spot in spots[i:i+savN-1]

      #Get thread ID
      id = Threads.threadid()

      #Copy mol to preserve original location
      new = deepcopy(mol)

      #Randomly rotate molecule (see hitAndStick for function code)
      randRotate!(new)

      #Spawn molecule and optimise 
      bdys = spawnMol(new, clu, com, spot, minD) |> (x -> opt(EoM, algo, x; pst...))

      #Get binding energy
      be = getPotEnergy(EoM, bdys) - (cluE + molE)

      #Save BE
      j = findfirst(e -> e == 0.0, ret[id, :])
      ret[id, j] += be
    end

    BE[i:i+savN-1] .+= vcat(ret...) |> dropZeros!
    df = mkBEdf(BE * -1)
    jldsave(inp["Saving"]["df"]; df)
    
    #increment i or change savN
    i + savN - 1 > N ? savN = N - i + 1 : i += savN

    #Rezero ret array
    ret .*= 0.0
  end

  df
end

function dropZeros!(x)
  i = findall(e -> e == 0.0, x)
  deleteat!(x, i)
end

function mkBEdf(arr)
  eVtoK  = 11604.525
  eVtocm = 8065.56

  df = Dict("eV" => arr, "K" => arr*eVtoK, "cm-1" => arr*eVtocm) |> DataFrame

  df
end

function getSpots(αShape)
  per = αShape.perimeter
  pts = αShape.pts
  
  ret = unique([j for i in per for j in i]) |> (x -> pts[x])

  [sum(pts[tri]) ./ 3 for tri in per] |> (x -> push!(ret, x...))

  return ret
end

function moveMol!(mol, v)
  for bdy in mol
    bdy.r .+= v
  end 
end

function spawnMol(mol, clu, com, spot, minD)

  vhat = (spot - com) ./ norm(spot - com)

  v    = spot + minD * vhat

  moveMol!(mol, v)
  
  d = [norm(j.r - i.r) for i in clu for j in mol] |> minimum

  while d < minD

    moveMol!(mol, vhat*0.1)

    d = [norm(j.r - i.r) for i in clu for j in mol] |> minimum

  end

  bdys = [clu; mol]

  return bdys
end