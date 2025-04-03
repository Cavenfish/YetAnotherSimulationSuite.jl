"""
Simple script to perform geometry optimization on a cell

Example TOML
-------------
[opt]
EoM = "MBX"
cell = "slab.xyz"
saveName = "slab_opt"
"""
using TOML, LinearAlgebra, JMD

inps = TOML.parsefile(ARGS[1])["opt"]

cell = readCell(inps["cell"])
EoM  = JMD.mkvar(inps["EoM"])

new = opt(EoM, JMD.ConjugateGradient(), cell, f_tol=1e-8, g_tol=1e-5, iterations=10000000)

writeCell("$(inps["saveName"]).xyz", new)

open("$(inps["saveName"]).txt", "w") do io
  eng  = getPotEnergy(EoM, new)
  fmax = getForces(EoM, new) |> (x -> norm.(x)) |> maximum
  println(io, "Energy: $(eng)")
  println(io, "Fmax  : $(fmax)")
end

