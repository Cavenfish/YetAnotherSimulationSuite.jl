"""
Calculates binding energies based on given
toml input card.

Example of Input Card
-----------------------

[Settings]
EoM = "COCO"
jldfile = "/home/brian/COjl/10K-MvH.jld2"
cluster = "500co"
sites = 250

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
  jld = cnfg["jldfile"]
  clu = cnfg["cluster"]
  N   = cnfg["sites"]

  # Load clusters
  jd = load(jld)

  # Pick cluster
  bdys = jd[clu]
end