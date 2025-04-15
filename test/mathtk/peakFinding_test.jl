
@testset "Peak and Turning Point Tests" begin
  f(x) = exp(- x^2 / 0.25)
  x    = collect(-2:0.005:2)
  y    = f.(x)
  i    = findPeaks(y)

  @test y[i][1] == 1.0

  g(x) = f(x+1.5) - 0.1*f(x+0.5) + f(x-1.5) - 0.1*f(x-0.5) + f(x)
  y    = g.(x)
  i    = findPeaks(y)

  @test x[i] == [-1.5, 0.0, 1.5]

  i    = findPeaks(y, min=0.95)

  @test x[i] == [-1.5, 1.5]

  i    = findPeaks(y, max=0.95)

  @test x[i][1] == 0.0

  i    = findTurningPoints(y)

  @test x[i] == [-1.5, -0.725, 0.0, 0.725, 1.5]
end