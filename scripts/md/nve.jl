"""
Simple script to run nve simulation

Example TOML
--------------
EoM = "HGNN"
bdys = "250co.xyz"
tmax = "50ps"
step = "fs"
saveName = "traj.jld2"
"""
using TOML, JLD2, JMD

inps = TOML.parsefile(ARGS[1])

EoM  = JMD.mkvar(inps["EoM"])
bdys = readXyz(inps["bdys"])
tmax = JMD.mkvar(inps["tmax"])
step = JMD.mkvar(inps["step"])

nve  = runMD(EoM, bdys, (0.0, tmax), step) |> processDynamics

jldsave(inps["saveName"]; traj=nve)