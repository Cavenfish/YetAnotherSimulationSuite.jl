"""
Anneals a cluster, based on given toml input card,
to get crystal structure.

Example of Input Card
-----------------------

[Settings]
EoM = "COCO"
jldfile = "/home/brian/COjl/10K-MvH.jld2"
cluster = "250co"
cycles = 10

[NVT]
thermostat = "BDP!"
thermoInps = "BDP"
temp = 25
kB = "kB"
tau = "100fs"
time = "5ps"

[OPT]
algo = "LBFGS()"

[Saving]
xyz = "crystal.xyz"
"""

function anneal(inpFile::String)
  inp = TOML.parsefile(inpFile)

  #Split dict for easier usage
  cnfg = inp["Settings"]
  nvtd = inp["NVT"]

  #Load up EoM and opt algo
  EoM  = mkvar(cnfg["EoM"])
  algo = mkvar(inp["OPT"]["algo"])

  #Prep nvt inputs
  T      = nvtd["temp"]
  tau    = mkvar(nvtd["tau"])
  nvtInp = mkvar(nvtd["thermoInps"])
  inps   = nvtInp(T, mkvar(nvtd["kB"]), tau)
  thermo = mkvar(nvtd["thermostat"])

  #Get leftover vars
  N = cfng["cycles"]
  t = mkvar(nvtd["time"])

  # Load clusters
  jd = load(jld)

  # Pick cluster
  bdys = jd[clu]

  #Start annealing
  for i in 1:N

    #Run NVT
    nvt = runNVT(EoM, (0, t), fs, bdys, thermo, inps)

    #Update bdys
    getLastFrame!(bdys, nvt)

    #Free memory
    @free nvt 

    #Optimize geometry
    bdys = opt(EoM, algo, bdys; iterations=100000)
    
  end

  writeXyz(inps["Saving"]["xyz"], bdys)

end