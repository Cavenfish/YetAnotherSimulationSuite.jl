@testset "LJ Calc Test" begin
  au = "../testingFiles/xyzFiles/Au20.xyz"

  d = Dict(
    "epsilon" => 0.2297,
    "sigma" => 2.95,
    "rs" => 19.0,
    "rc" => 20.0
  )

  file = joinpath(@__DIR__, au)
  bdys = readSystem(file)
  calc = LJ(d)
  new  = opt(calc, ConjugateGradient(), bdys; g_tol=1e-9)

  E = getPotEnergy(calc, new)
  F = getForces(calc, new)
  f = norm.(F) |> maximum

  @test isapprox(E, -15.541123; atol=1e-5)
  @test f < 1e-8
end
