
@testset "Test Geometric Operations" begin
  f    = "../testingFiles/xyzFiles/co.xyz"
  file = joinpath(@__DIR__, f)
  co   = readSystem(file)

  translate!(co, [0.0, 1.0, 0.0])

  # check translation worked
  @test co[1].r == [0.0, 1.0, 0.0]
  @test co[2].r == [1.1, 1.0, 0.0]

  old = CoM(co)
  randRotate!(co, old)

  # check CoM is unchanged by rotation about it
  @test isapprox(CoM(co), old; atol=1e-8)

  co = readSystem(file)

  # rotate 90 degrees around y-axis
  rotate!(co, (0.0, -pi/2, 0.0))

  @test isapprox(co[1].r, [0.0,0.0,0.0]; atol=1e-10)
  @test isapprox(co[2].r, [0.0,0.0,1.1]; atol=1e-10)

  # rotate 90 degrees around y-axis about [0,0,0.55]
  rotate!(co, (0.0, pi/2, 0.0), [0.0,0.0,0.55])

  @test isapprox(co[1].r, [-0.55,0.0,0.55]; atol=1e-10)
  @test isapprox(co[2].r, [0.55,0.0,0.55]; atol=1e-10)

  # Move CoM to origin then rotate about origin
  old = CoM(co)
  translate!(co, -old)
  randRotate!(co)

  # check CoM is unchanged by rotation
  @test isapprox(CoM(co), [0.0,0.0,0.0]; atol=1e-8)
end