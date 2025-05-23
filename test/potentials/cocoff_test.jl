
@testset "Dimer Test" begin

  # Load dimer
  bdys = joinpath(@__DIR__, "../testingFiles/xyzFiles/co-dimer.xyz") |> readSystem

  # Prep calculator
  calc = MvHff()

  # Optimize dimer
  new = opt(calc, JMD.LBFGS(), bdys, iterations=1000, g_tol=1e-8)

  # Get distances and interaction energy
  r1 = new[1].r - new[2].r |> norm
  r2 = new[3].r - new[4].r |> norm
  r3 = CoM(new[1:2]) - CoM(new[3:4]) |> norm
  E  = getPotEnergy(calc, new) / 0.000124

  # Test statements
  @test r1 ≈ 1.128 atol=5e-4
  @test r2 ≈ 1.128 atol=5e-4
  @test r3 ≈ 3.7   atol=5e-2
  @test E  ≈ -128  atol=5e-1
end
