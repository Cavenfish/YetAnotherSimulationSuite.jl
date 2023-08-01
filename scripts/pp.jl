md = include("./src/JMD.jl")

jd = md.load("./10K-MvH.jld2")

bdys = jd["50co"]

nve = md.runNVE(md.COCO, (0, 1*md.ps), md.fs, bdys)

@time traj = md.processDynamics(nve)

@time traj = md.processDynamics(nve)

@time traj = md.processDynamics(nve)
