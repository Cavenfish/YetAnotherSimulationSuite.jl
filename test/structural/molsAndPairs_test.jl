
@testset "Distance Matrix" begin
  bdys::Vector{MyAtoms} = [
    Particle([1.000,  0.000, 0.000], zeros(3), 15.999, "O"),
    Particle([1.500,  0.000, 0.000], zeros(3), 15.999, "O")
  ]

  lat  = [2.0 0 0; 0 2.0 0; 0 0 2.0]
  cell = makeCell(bdys, lat)

  D1 = YASS.distanceMatrixOrthorhombicCell(cell)
  D2 = YASS.distanceMatrixAnyCell(cell)

  @test D1[1, 2] == 0.5
  @test D2[1, 2] == 0.5
  @test D1[2, 1] == 0.5
  @test D2[2, 1] == 0.5
end