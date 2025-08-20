
@testset "TEST VACF" begin
  file = joinpath(@__DIR__, "../testingFiles/xyzFiles/h2o.xyz")
  bdys = readSystem(file)
  calc = TIP4Pf()
  new  = opt(calc, YASS.LBFGS(), bdys; g_tol=1e-10, iterations=500)
  v,m  = getHarmonicFreqs(calc, new)

  vibExcite!(new, m[:, end], 0.05)

  tj       = run(TIP4Pf(), new, (0.0, 15ps), 0.1fs, YASS.NVE())
  vel, mas = getVelMas(tj)
  inps     = vacfInps(vel, mas, 1e16, true, YASS.HannM, 8, true)
  out      = VDOS(inps)

  # Get peaks to test they match with water peak
  m = maximum(out.I) * 0.95
  i = findPeaks(out.I; min=m)[1]
  D = getDiffusionCoefficient(out)

  @test isapprox(out.v[i], 3687; atol=5)
  @test isapprox(D, 0.0; atol=1)
end
