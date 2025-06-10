
@testset "Test Particles" begin
  f    = "../testingFiles/xyzFiles/h2o.xyz"
  file = joinpath(@__DIR__, f)
  h2o  = readSystem(file)

  centerBdys!(h2o)

  @test isapprox(CoM(h2o), [0.0,0.0,0.0]; atol=1e-6)

  a = h2o[2].r |> deepcopy
  b = h2o[3].r |> deepcopy

  translateBdys!(h2o, [1,0,0])
  swapAtoms!(h2o, 2, 3)

  @test isapprox(a, h2o[3].r .- [1,0,0]; atol=1e-8)
  @test isapprox(b, h2o[2].r .- [1,0,0]; atol=1e-8)

  swapIso!(h2o, [2,3], [2.0, 2.0])

  @test h2o[2].m == 2.0
  @test h2o[3].m == 2.0

end
