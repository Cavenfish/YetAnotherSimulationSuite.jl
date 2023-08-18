
function vibDisp(inpFile::String)

  # Read in input card
  f   = open(inpFile, "r")
  inp = TOML.parse(f)
  close(f)

  # Split dict for easier usage
  expt = inp["expt"]
  savi = inp["saving"]

  # Load up EoM 
  fn  = Symbol(expt["EoM"])
  EoM = @eval $fn

  # Load other vars for easier usage
  E      = expt["energy"]
  time   = expt["time"] * ps
  loc    = expt["location"]
  iso    = expt["isotope"]
  jld    = expt["jldfile"]
  clu    = expt["cluster"]
  splits = expt["splits"]

  # Load clusters
  jd = load(jld)

  # Pick cluster
  bdys = jd[clu]
  zeroVCoM!(bdys)

  #Swap cluster isotope
  if "cluIso" in keys(expt)
    tmp  = expt["cluIso"]
    swap = collect(1:length(bdys))
    mas  = [tmp[1], tmp[2] for i in 1:length(bdys)/2]
    swapIso!(bdys, swap, mas)
  end

  # Randomly select molecule to excite
  mol = pickRandomMol(bdys, loc)

  # Swap mass of CO
  swapIso!(bdys, mol, iso)

  # Run short NVE to equilibrate system
  equil = runNVE(EoM, (0, 10ps), fs, bdys)
  getLastFrame!(bdys, equil)

  # Save equil data
  open("0.tmp", "w") do f
    serialize(f, equil)
  end

  # Excite CO
  co  = bdys[mol]
  f,m = getHarmonicFreqs(EoM, co)
  vibExcite!(co, m[:,6], E)

  # Run in parts to avoid segfaults
  for i in 1:splits
    t1  = ((i-1)/splits) * time + 10ps
    t2  = (i/splits) * time + 10ps
    nve = runNVE(EoM, t1, t2, fs, bdys)
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
  df   = trackEnergyDissipation(traj, EoM, mol)

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

end