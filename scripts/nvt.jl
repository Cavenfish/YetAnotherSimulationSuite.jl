include("./src/JMD.jl")

bdys = JMD.readXyz("./xyzFiles/100.xyz")

thermoInps = JMD.BDP(10.0, JMD.kB, 100*JMD.fs)

@time solu = JMD.runNVT(JMD.HGNN, (0, 5*JMD.fs), JMD.fs, bdys, JMD.BDP!, thermoInps)

@time solu = JMD.runNVT(JMD.HGNN, (0, 5000*JMD.fs), JMD.fs, bdys, JMD.BDP!, thermoInps)

#JMD.writeXyzTraj("20coNVT.xyz", solu)
