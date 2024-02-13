"""
Runs a vibrational energy dissipation simulation 
based on conditions given via a toml input card.

Example of Input Card 
-----------------------

[expt]
EoM = "COCO"
energy = 0.4
time = 2000
dt = 1
location = "bulk"
isotope = [12.011, 15.999]
jldfile = "/home/brian/COjl/10K-MvH.jld2"
cluster = "250co"
splits = 100
cluIso = [13.003, 17.999]
KE = 10
vzpe.energy = 0.001
vzpe.modes = [1,2,3]

[vacf]
inter = 10
safe = 15000

[track]
Rs = [[3,5],[5,9],[9,30]]

[saving]
df = "co-am13C18O_MvH_DF_1.jld2"
tj = "co-am13C18O_MvH_TJ_1.jld2"
vd = "co-am13C18O_MvH_VD_1.jld2"
re = "co-am13C18O_MvH_RE_1.jld2"
"""

function vibDisp(inpFile::String)

  # Read in input card
  f   = open(inpFile, "r")
  inp = TOML.parse(f)
  close(f)

  # Split dict for easier usage
  expt = inp["expt"]
  savi = inp["saving"]

  # Load timestep
  "dt" in keys(expt) ? dt = expt["dt"]*fs : dt = fs

  # Load up EoM 
  fn  = Symbol(expt["EoM"])
  EoM = @eval $fn

  # Load other vars for easier usage
  E      = expt["energy"]
  time   = expt["time"] * ps
  loc    = expt["location"]
  iso    = expt["isotope"]
  jld    = expt["jldfile"]
  splits = expt["splits"]

  # Load clusters
  jd = load(jld)

  # Pick cluster
  bdys = if "cluster" in keys(expt)
           clu = expt["cluster"]       
           jd[clu]
         else
           k = keys(jd)
           jd[rand(k)]
         end
  zeroVCoM!(bdys)

  #Swap cluster isotope
  if "cluIso" in keys(expt)
    tmp  = expt["cluIso"]
    swap = collect(1:length(bdys))
    mas  = repeat(tmp, div(length(bdys),2) )
    swapIso!(bdys, swap, mas)
  end

  # Randomly select molecule to excite
  mol = pickRandomMol(bdys, loc) |> (x -> findall(e -> e in x, bdys))

  # Swap mass of CO
  swapIso!(bdys, mol, iso)

  #Add VZPE to cluster
  if "vzpe" in keys(expt)
    zpe = expt["vzpe"]["energy"]
    f,m = getHarmonicFreqs(EoM, bdys)
    N   = div(length(bdys), 2)
    m   = m[:, end-(N-1):end]
    for i = expt["vzpe"]["modes"]
      vibExcite!(bdys, m[:, i], zpe)
    end
    expt["vzpe"]["modes"] = [m[:, i] for i in expt["vzpe"]["modes"]]
  end

  # Run short NVE to equilibrate system
  equil = runNVE(EoM, (0, 10ps), dt, bdys, save="sparse")
  getLastFrame!(bdys, equil)

  # Save equil data
  open("0.tmp", "w") do f
    serialize(f, equil)
  end

  # Excite CO
  co  = bdys[mol]
  f,m = getHarmonicFreqs(EoM, co)
  vibExcite!(co, m[:,6], E)

  #Translational excitation
  if "KE" in keys(expt)
    KE = expt["KE"]
    transExcite!(co, KE)
  end

  # Run in parts to avoid segfaults
  for i in 1:splits
    t1  = ((i-1)/splits) * time + 10ps
    t2  = (i/splits) * time + 10ps
    nve = runNVE(EoM, t1, t2, dt, bdys, save="sparse")
    getLastFrame!(bdys, nve)
    
    open("$i.tmp", "w") do f
      serialize(f, nve)
    end

    #Free memory
    @free nve
  end

  # Post-process
  tmp  = ["$i.tmp" for i in 0:splits]
  traj = processTmpFiles(tmp; step=100)
  df   = if "vzpe" in keys(expt)
    trackEnergyDissipation(traj, EoM, mol, eignvec=expt["vzpe"]["modes"])
  else
    trackEnergyDissipation(traj, EoM, mol)
  end

  # Load saving vars for easy usage
  dfName = savi["df"] 
  tjName = savi["tj"]

  # Save data
  jldsave(dfName; df)
  jldsave(tjName; traj)

  if "vacf" in keys(inp)
    vdName = savi["vd"]
    safe   = inp["vacf"]["safe"]

    j   = inp["vacf"]["inter"]
    tmp = ["$i.tmp" for i in 1:j:splits]

    vd  = trackVACF(tmp, safe)

    jldsave(vdName; vd)
  end

  if "track" in keys(inp)
    reName = savi["re"]
    Rs     = inp["track"]["Rs"]

    re = trackRadialEnergy(traj; pot=EoM, Rs=Rs)

    jldsave(reName; re)
  end

end