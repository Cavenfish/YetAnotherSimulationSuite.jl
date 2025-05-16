
@testset "Test Harmonic Frequency" begin
  file = joinpath(@__DIR__, "../testingFiles/xyzFiles/co-dimer.xyz")
  bdys = readSystem(file)
  
  # ensure bdys are optimized
  new = opt(MvHff, JMD.LBFGS(), bdys; g_tol=1e-9, iterations=500)

  v,m = getHarmonicFreqs(MvHff, new)

  @test length(v) == 3*length(bdys)
  @test size(m) == (length(v), length(v))
  @test isapprox(real(v[end]), 2197, atol=1)

  animateMode(bdys, m[:, end], "./tmp.xyz")

  # not a great test but the above execution at least
  # lets me know code changes break the function call
  @test isfile("./tmp.xyz")
  rm("./tmp.xyz")

  @testset "Test IPR and PR" begin
    pr = getIPR(m) |> (y -> getPR(y; x=0.25))

    @test pr[end] == 1
  end
end