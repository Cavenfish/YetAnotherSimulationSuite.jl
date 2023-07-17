
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
  E    = expt["energy"]
  time = expt["time"] * ps
  loc  = expt["location"]
  iso  = expt["isotope"]
  jld  = expt["jldfile"]
  clu  = expt["cluster"]

  # Load clusters
  jd = load(jld)

  # Pick cluster
  bdys = jd[clu]
  zeroVCoM!(bdys)

  # Randomly select molecule to excite
  mol = pickRandomMol(bdys, loc)

  # Swap mass of CO
  swapIso!(bdys, mol, iso)

  # Run short NVE to equilibrate system
  equil = runNVE(EoM, (0, 10ps), fs, bdys)
  getLastFrame!(bdys, equil)

  # Free Memory
  equil = 0
  GC.gc()

  # Excite CO
  co  = bdys[mol]
  f,m = getHarmonicFreqs(EoM, co)
  vibExcite!(co, m[:,6], E)

  # Run NVE
  nve  = runNVE(EoM, (0, time), fs, bdys)

  # Post-process 
  traj = processDynamics(nve; step=100)
  df   = trackEnergyDissipation(traj, EoM, mol)
  v,m  = getVelMas(nve)

  # For now, wont be kept
  writeXyzTraj(tjName * ".xyz", nve; dt=100)

  # Free Memory
  nve = 0
  GC.gc()

  # Load saving vars for easy usage
  dfName = savi["df"] 
  tjName = savi["tj"]
  vmName = savi["vm"]

  # Save data
  jldsave(dfName; df)
  jldsave(tjName * ".jld2"; traj)
  jldsave("./myVnM.jld2"; v, m)

end