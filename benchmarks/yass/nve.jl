using Pkg
Pkg.activate("./")

using YetAnotherSimulationSuite

i = ARGS[1]
j = parse(Float64, ARGS[2])

bdys = readSystem("../xyzfiles/$(i)au.xyz")

d = Dict(
  "epsilon" => 0.2297,
  "sigma" => 2.95,
  "rs" => 19.0,
  "rc" => 20.0
)

calc = LJ(d)

nve = run(calc, bdys, j * 1u"ps", 1u"fs", NVE())