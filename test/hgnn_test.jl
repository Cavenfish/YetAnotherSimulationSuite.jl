using Test
using LinearAlgebra

include("../src/JMD.jl")

@testset "HGNN Test" begin

  o1   =  [ 0.507464,   -5.559636,    0.000000] ./ 0.5291772083
  c1   =  [-0.676403,   -4.827579,    0.000000] ./ 0.5291772083
  o2   =  [-0.002646,    0.507559,   -0.026816] ./ 0.5291772083
  c2   =  [ 0.003527,   -0.676530,    0.035743] ./ 0.5291772083
  
  r     = [c1,o1,c2,o2]
  F     = [zeros(Float64, 3) for i in 1:4]
  rhats = zeros(Float64, 6, 3)
  hats  = zeros(Float64, 6, 3)
  dPdr  = zeros(Float64, 7, 6)
  P     = zeros(Float64, 7)
  A     = zeros(Float64, 7, 45)
  dv    = zeros(Float64, 6)

  JMD.getUnitVectors!(hats, r[1],r[2],r[3],r[4])

  vco1   = JMD.molPot!(F, r, [1,2], JMD.hgnnMolVars)
  vco2   = JMD.molPot!(F, r, [3,4], JMD.hgnnMolVars)

  # pull magnitude of force
  diff   = r[2] .- r[1]
  R      = sqrt(dot(diff, diff))
  v      = diff / R
  dvco1  = dot(F[1], v)

  # pull magnitude of force
  diff   = r[4] .- r[3]
  R      = sqrt(dot(diff, diff))
  v      = diff / R
  dvco2  = dot(F[3], v)

  # pull magnitude of forces
  vint   = JMD.pairPot!(F, r, [1,2] => [3,4], JMD.hgnnPairVars, rhats, dPdr, P, A)
  for i in 1:6
    dv[i] = dot(rhats[i,:], hats[i,:])
  end

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
