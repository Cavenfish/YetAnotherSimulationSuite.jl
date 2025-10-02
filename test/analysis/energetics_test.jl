
@testset "Test Energetics Functions" begin
  file = joinpath(@__DIR__, "../testingFiles/xyzFiles/h2o.xyz")
  bdys = readSystem(file)
  calc = SPC("SPC/F")
  new1 = opt(calc, LBFGS(), bdys; g_tol=1e-9, iterations=500)
  new2 = deepcopy(new1)

  @test isapprox(0.0, getRotEnergy(new1); atol=1e-8)
  @test isapprox(0.0, getTransEnergy(new1); atol=1e-8)
  transExcite!(new1, 0.1)
  @test isapprox(0.0, getRotEnergy(new1); atol=1e-8)
  @test isapprox(0.1, getTransEnergy(new1); atol=1e-5)

  v,m  = getHarmonicFreqs(calc, new2)

  @test isapprox(0.0, getVibEnergy(new2, m[:, end]; calc=calc); atol=1e-10)
  vibExcite!(new2, m[:, end], 0.1)
  @test isapprox(0.0, getTransEnergy(new2); atol=1e-8)
  @test isapprox(0.1, getVibEnergy(new2, m[:, end]; calc=calc); atol=1e-5)

  vcom = vCoM(new2) |> norm
  @test isapprox(0.0, vcom; atol=1e-8)

end
