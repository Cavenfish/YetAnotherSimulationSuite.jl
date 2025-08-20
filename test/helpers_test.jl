
@testset "CoM and vCoM Tests" begin
  μ    = ((1/5) + (1/10))^-1
  com  = [35/15, 0.0, 0.0]
  vcom = [0.0, 5/15, 10.0/15]
  bdys = [
    Particle([1.0,0.0,0.0], [0.0,1.0,0.0], 5.0, "X"),
    Particle([3.0,0.0,0.0], [0.0,0.0,1.0], 10.0, "Y")
  ]
  pos  = [i.r for i in bdys]
  vel  = [i.v for i in bdys]
  mas  = [i.m for i in bdys]

  @test YASS.reducedMass(bdys) == μ
  @test YASS.reducedMass(mas) == μ
  @test CoM(bdys ) == com
  @test vCoM(bdys) == vcom
  @test CoM(pos,mas) == com
  @test vCoM(vel,mas) == vcom

  zeroVCoM!(bdys)

  @test vCoM(bdys) == [0.0, 0.0, 0.0]

end