#md = include("./src/JMD.jl")
md = include("../../Code/JMD/src/JMD.jl")

mol = "rotTest/1.xyz"

thermoInp = md.BDP(10.0, md.kB, 0.1*md.ps)

inp = md.HnS(20, 5, 5, 10.0, mol, mol, "tmp.xyz", md.BDP!, thermoInp)

@time md.hitAndStick(md.HGNNdyn, inp)

#inp = md.HnS(20, 5, 5, 10.0, mol, mol, "tmp.xyz", md.BDP!, thermoInp)

@time md.hitAndStick(md.COCOdyn, inp)

#inp = md.HnS(40, 5, 5, 10.0, mol, mol, "tmp.xyz", md.BDP!, thermoInp)

#@time md.hitAndStick(md.COCOdyn, inp)

#inp = md.HnS(80, 5, 5, 10.0, mol, mol, "tmp.xyz", md.BDP!, thermoInp)

#@time md.hitAndStick(md.COCOdyn, inp)
