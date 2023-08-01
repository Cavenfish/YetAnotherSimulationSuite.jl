include("./src/JMD.jl")

dt   = JMD.fs # 1fs in ASE time

#bdys = JMD.readASExyz("../CO_Project/xyz_files/20co.xyz")
bdys = JMD.readXyz("./xyzFiles/10_opt.xyz")

@time solu = JMD.runNVE(JMD.COCO, (0, 1*dt), dt, bdys)

@time solu = JMD.runNVE(JMD.COCO, (0, 2000*dt), dt, bdys)

@time solu = JMD.runNVE(JMD.COCO, (0, 2000*dt), dt, bdys)

#@time solu = JMD.runNVE(JMD.COCOtst, (0, 1*dt), dt, bdys)

#@time solu = JMD.runNVE(JMD.COCOtst, (0, 2000*dt), dt, bdys)

#@time solu = JMD.runNVE(JMD.COCOtst, (0, 2000*dt), dt, bdys)
