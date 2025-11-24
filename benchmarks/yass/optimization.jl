using Pkg
Pkg.activate("./")

using Optim
using YetAnotherSimulationSuite

i = ARGS[1]
j = parse(Int, ARGS[2])

bdys = readSystem("../xyzfiles/$(i)au.xyz")

d = Dict(
  "epsilon" => 0.2297,
  "sigma" => 2.95,
  "rs" => 19.0,
  "rc" => 20.0
)

calc = LJ(d)

new = opt(calc, LBFGS(), bdys; g_tol=1e-20, iterations=j)

