using Pkg, TOML
Pkg.activate("../")
md = include("../src/JMD.jl")

#Read in inputs
f   = open("./disp.toml", "r")
inp = TOML.parse(f)
close(f)

# Predefine some things (ARGS eventually)
EoM  = md.COCO
E    = 0.4
time = 20*md.ps
loc  = "bulk"
iso  = [12.011, 17.999]

# Load clusters
jd = md.load("../10K-MvH.jld2")

# Pick cluster
bdys = jd["500co"]
md.zeroVCoM!(bdys)

# Randomly select molecule to excite
mol = md.pickRandomMol(bdys, loc)

# Swap mass of CO
md.swapIso!(bdys, mol, iso)

# Run short NVE to equilibrate system
println("Equil System")
equil = md.runNVE(EoM, (0, 10*md.ps), md.fs, bdys)
md.getLastFrame!(bdys, equil)

# Free Memory
equil = 0
GC.gc()

# Excite CO
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
v,m  = md.getVelMas(nve)

# Save data
println("Save Data")
md.jldsave("./myData.jld2"; df)
md.jldsave("./myTraj.jld2"; traj)
md.jldsave("./myVnM.jld2"; v, m)

# For now, wont be kept
md.writeXyzTraj("./myTraj.xyz", nve; dt=100)
