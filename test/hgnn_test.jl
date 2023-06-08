using Test

include("../src/JMD.jl")

@testset "HGNN Test" begin

  o1   =  [ 0.507464,   -5.559636,    0.000000] ./ 0.5291772083
  c1   =  [-0.676403,   -4.827579,    0.000000] ./ 0.5291772083
  o2   =  [-0.002646,    0.507559,   -0.026816] ./ 0.5291772083
  c2   =  [ 0.003527,   -0.676530,    0.035743] ./ 0.5291772083
  # vars = JMD.readInVars("../src/md/potentials/params/nn_ococ_w20.txt")
  dPdr  = zeros(Float64, 7, 6)
  P     = zeros(Float64, 7)

  vint, dv = JMD.pairPot([c1,o1], [c2,o2], JMD.hgnnPairVars, dPdr, P)
  vco1, dvco1, f = JMD.molPot([c1,o1], JMD.hgnnMolVars)
  vco2, dvco2, f = JMD.molPot([c2,o2], JMD.hgnnMolVars)

  vtot = vint + vco1 + vco2
  dv[1] += dvco1
  dv[2] += dvco2

  @test vint ≈ -76.2   atol=1e-1
  @test vtot ≈ 19711.9 atol=1e-1

  @test dv[1] ≈ 53917.4 atol=1e-1
  @test dv[2] ≈ 23683.8 atol=1e-1
  @test dv[3] ≈ 40.2 atol=1e-1
  @test dv[4] ≈ -26.8 atol=1e-1
  @test dv[5] ≈ -21.9 atol=1e-1
  @test dv[6] ≈ 50.9 atol=1e-1

end
