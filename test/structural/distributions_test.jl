
@testset "RDF Tests" begin
  file = joinpath(@__DIR__, "../testingFiles/xyzFiles/iceIh.xyz")
  bdys = readSystem(file) |> makeBdys

  # Molecule CoM RDF
  x, y = rdf(bdys)
  Ro   = x[findPeaks(y, min=0.04)]

  @test isapprox(Ro[1], 2.70, atol=5e-2)
  @test isapprox(Ro[2], 4.45, atol=5e-2)

  # Oxygen atom RDF
  x, y = rdf(bdys, "O")
  Ro   = x[findPeaks(y, min=0.04)]

  @test isapprox(Ro[1], 2.70, atol=5e-2)
  @test isapprox(Ro[2], 4.45, atol=5e-2)

  # O-H RDF
  x, y = rdf(bdys, "O" => "H")
  Ro   = x[findPeaks(y, min=0.04)]

  @test isapprox(Ro[1], 0.91, atol=5e-2)
end