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

  #Get αShape of cluster
  αShape = [CoM(bdys[i:i+1]) for i in 1:2:length(bdys)] |> alphashape

  #Get BE for N number of sites
  for spot in αShape[rand(1:length(αShape), N)]

    #Randomly rotate molecule (see hitAndStick for function code)
    randRotate!(mol)

    #Spawn molecule and optimise 
    opt = spawnMol(mol, clu, αShape, minD) |> (x -> opt(EoM, algo, x; kwargs...))

    #Get binding energy
    be = getPotEnergy(EoM, opt) - (cluE + molE)

  end

end

function spawnMol(mol, clu, spot, minD)

  return bdys
end