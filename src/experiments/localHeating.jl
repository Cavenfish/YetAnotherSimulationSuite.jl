"""
Simulates local heating based on given toml input card.

Example of Input Card
-----------------------

[expt]
EoM = "COCO"
clu = "/home/brian/COjl/clu.jld2"
cluName = "bdys"
time = 500

[NVT]
tau = 100
time = 100
radius = 5
temp = [10, 75]

[saving]
tj = "result.jld2"
"""

function localHeat(inpFile::String)
  inp = TOML.parsefile(inpFile)

  #Split dict for easier usage
  expt = inp["expt"]
  nvtI = inp["NVT"]

  #Load up EoM
  EoM = mkvar(expt["EoM"])

  #Load cluster
  bdys = jldopen(expt["clu"])[expt["cluName"]]

  #Randomly pick central atom
  mol = rand(bdys)

  #Find all atoms to be locally heated
  R2 = findall(e -> norm(e.r - mol.r) <= 5, bdys)
  R1 = [i for i=1:length(bdys) if !(i in R2)]

  #Prep NVT inps
  thermoInps = BDPnT(nvtI["temp"], [R1, R2], kB, nvtI["tau"] * fs)

  #Run two-temp NVT
  nvt = md.runNVT(EoM, (0, nvtI["time"]*ps), fs, bdys, BDPnT!, thermoInps; save="sparse")

  #Get Last frame and zero out vCoM
  new = getLastFrame(nvt) |> zeroVCoM

  #Run NVE
  nve = runNVE(EoM, (0, expt["time"] * ps), fs, new; save="sparse") 

  #Process dynamics
  tj = processDynamics(nve)

  #Save traj
  jldsave(inp["saving"]["tj"]; tj)
end