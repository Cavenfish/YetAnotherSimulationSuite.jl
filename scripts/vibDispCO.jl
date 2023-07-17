using Pkg
Pkg.activate("../")
md = include("../src/JMD.jl")

# Predefine some things (ARGS eventually)
EoM  = md.HGNN
E    = 0.4
time = 2*md.ps
loc  = "bulk"
iso  = [12.011, 17.999]

# Load clusters
jd = md.load("../10K-HGNN.jld2")

# Pick cluster
bdys = jd["50co"]
md.zeroVCoM!(bdys)

# Randomly select molecule to excite
mol = md.pickRandomMol(bdys, loc)

# Swap mass of CO
md.swapIso!(bdys, mol, iso)

# Run short NVE to equilibrate system
println("Equil System")
equil = md.runNVE(EoM, (0, 1*md.ps), md.fs, bdys)
md.getLastFrame!(bdys, equil)

# # Excite CO
println("Excite CO")
co  = bdys[mol]
f,m = md.getHarmonicFreqs(EoM, co)
md.vibExcite!(co, m[:,6], E)

# Run NVE
println("Run NVE")
nve  = md.runNVE(EoM, (0, time), md.fs, bdys)

println("Post-process")
traj = md.processDynamics(nve; step=100)
df   = md.trackEnergyDissipation(traj, EoM, mol)

# Save data
println("Save Data")
md.jldsave("./myData3.jld2"; df)
md.jldsave("./myTraj3.jld2"; traj)

# For now, wont be kept
md.writeXyzTraj("./myTraj3.xyz", nve; dt=100)