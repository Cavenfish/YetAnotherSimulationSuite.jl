using Pkg
Pkg.activate("./")

using YetAnotherSimulationSuite

i = ARGS[1]

bdys = readSystem("../xyzfiles/$(i)au_opt.xyz")

d = Dict(
  "epsilon" => 0.2297,
  "sigma" => 2.95,
  "rs" => 19.0,
  "rc" => 20.0
)

calc = LJ(d)

v,m = getHarmonicFreqs(calc, bdys)