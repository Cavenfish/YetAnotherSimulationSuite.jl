
@testset "Stress Tensor" begin
  fnam = joinpath(@__DIR__, "../testingFiles/xyzFiles/iceIh_small.xyz")
  cell = readSystem(fnam)

  g1 = getNumericalStress(SPCF(), cell)
  g2 = getNumericalStressOrthogonal(SPCF(), cell)

  @test g1[1] == g2[1]
  @test g1[5] == g2[5]
  @test g1[9] == g2[9]

  @test isapprox(g1[1], -0.0420102; atol=1e-7)
  @test isapprox(g1[2], -0.000687213; atol=1e-7)
  @test isapprox(g1[3], -9.48539e-5; atol=1e-7)
end