"""
Calculates binding energies based on given
toml input card.

Example of Input Card
-----------------------

[Settings]
EoM = "COCO"
mol = "/home/brian/COjl/xyzFiles/1.xyz"
jldfile = "/home/brian/COjl/10K-MvH.jld2"
cluster = "500co"
sites = 250
minD = 3.5

[OPT]
algo = "LBFGS"
kwargs.x_tol = 1e-4
kwargs.f_tol = 1e-12
kwargs.g_tol = 1e-4
kwargs.iterations = 100000

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

  #Build opt kwargs dict
  kwargs = Dict()
  for key in keys(inp["OPT"]["kwargs"])
    value = inp["OPT"]["kwargs"][key]
    kwargs[Symbol(key)] = value
  end

  #Get leftover vars
  minD = cnfg["minD"]
  jld  = cnfg["jldfile"]
  N    = cnfg["sites"]

  # Load clusters
  jd = load(jld)

  #Read in and pre-optimise cluster and mol
  mol = readXyz(cnfg["mol"]) |> (x -> opt(EoM, algo, x; kwargs...))
  clu = jd[cnfg["cluster"]]  |> (x -> opt(EoM, algo, x; kwargs...))

  #Get cluster center of mass (for surface norms later)
  com = CoM(clu)

  # Get energies
  molE = getPotEnergy(EoM, mol)
  cluE = getPotEnergy(EoM, clu)

  #Get spots on surface of cluster
  spots = [CoM(clu[i:i+1]) for i in 1:2:length(clu)] |> alphashape |> getSpots

  #Check N ≤ spots
  N ≤ length(spots) || throw("Not enough spots found for $N sites")

  #Prep BE df
  ret = Float64[]

  #Get BE for N number of sites
  for spot in sample(spots, N; replace=false)

    #Copy mol to preserve original location
    new = deepcopy(mol)

    #Randomly rotate molecule (see hitAndStick for function code)
    randRotate!(new)

    #Spawn molecule and optimise 
    bdys = spawnMol(new, clu, com, spot, minD) |> (x -> opt(EoM, algo, x; kwargs...))

    #Get binding energy
    be = getPotEnergy(EoM, bdys) - (cluE + molE)

    push!(ret, be)
  end

  jldsave(inp["Saving"]["df"]; ret)

  return ret
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