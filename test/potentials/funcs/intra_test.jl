rbuf = zeros(3)
Fbuf = zeros(3)

@testset "Morse Test" begin
  u = [[0.0,0.0,0.0], [1.0,0.0,0.0]]
  F = zero(u)
  Z = zeros(3)
  
  # Test equilibrium
  E1    = JMD._Morse!(F, u, 1, 2, 50.0, 2.2, 1.0)
  E2, f = JMD._Morse(1.0, u[2], 50.0, 2.2, 1.0)
  @test E1 == 0.0
  @test E2 == 0.0
  @test F[1] == Z
  @test F[2] == Z
  @test f == Z
  
  E1    = JMD._Morse!(F, u, Fbuf, rbuf, 1, 2, 50.0, 2.2, 1.0)
  @test E1 == 0.0
  @test F[1] == Z
  @test F[2] == Z

  # Test non-equilibrium
  E  = 15.0 * (1 - exp(-2.0 * (2.0 - 1.0)))^2
  Z .= -30.0 * 2.0 * exp(-2.0 * (2.0 - 1.0)) * (1 - exp(-2.0 * (2.0 - 1.0))) .* [1.0, 0.0, 0.0] 
  
  u[2] .+= [1.0, 0.0, 0.0]
  E1    = JMD._Morse!(F, u, 1, 2, 15.0, 2.0, 1.0)
  E2, f = JMD._Morse(2.0, u[2], 15.0, 2.0, 1.0)
  @test E1 == E
  @test E2 == E
  @test f == Z
  @test F[1] == -Z
  @test F[2] == Z

  F     = zero(u)
  E1    = JMD._Morse!(F, u, Fbuf, rbuf, 1, 2, 15.0, 2.0, 1.0)
  @test E1 == E
  @test F[1] == -Z
  @test F[2] == Z

end

@testset "Harmonic Bond Test" begin
  u = [[0.0,0.0,0.0], [1.0,0.0,0.0]]
  F = zero(u)
  Z = zeros(3)

  # Test equilibrium
  E1    = JMD._harmonicBond!(F, u, 1, 2, 50.0, 1.0)
  E2, f = JMD._harmonicBond(1.0, u[2], 50.0, 1.0)
  @test E1 == 0.0
  @test E2 == 0.0
  @test F[1] == Z
  @test F[2] == Z
  @test f == Z

  E1    = JMD._harmonicBond!(F, u, Fbuf, rbuf, 1, 2, 50.0, 1.0)
  @test E1 == 0.0
  @test F[1] == Z
  @test F[2] == Z

  # Test non-equilibrium
  E  = 0.5 * 50.0 * (2.0 - 1.0)^2 
  Z .= - 50.0 * (2.0 - 1.0) .* [1.0,0.0,0.0]

  u[2] .+= [1.0, 0.0, 0.0]
  E1    = JMD._harmonicBond!(F, u, 1, 2, 50.0, 1.0)
  E2, f = JMD._harmonicBond(2.0, u[2], 50.0, 1.0)
  @test E1 == E
  @test E2 == E
  @test f == Z
  @test F[1] == -Z
  @test F[2] == Z

  F     = zero(u)
  E1    = JMD._harmonicBond!(F, u, Fbuf, rbuf, 1, 2, 50.0, 1.0)
  @test E1 == E
  @test F[1] == -Z
  @test F[2] == Z
end

@testset "Harmonic Angle Test" begin
  u = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0]]
  F = zero(u)
  Z = zero(u)

  # Test equilibrium
  E1             = JMD._harmonicBondAngle!(F, u, 2, 1, 3, 50.0, pi/2)
  E2, f1, f2, fo = JMD._harmonicBondAngle(u[2]-u[1], u[3]-u[1], 50.0, pi/2)
  @test E1 == 0.0
  @test E2 == 0.0
  @test F[1] == Z[1]
  @test F[2] == Z[2]
  @test F[3] == Z[3]
  @test f1 == f2 == fo == Z[1]

  # Test non-equilibrium
  u[3] .= [cos(3pi/4),sin(3pi/4),0.0]
  E     = 0.5 * 50.0 * (3pi/4 - pi/2)^2
  ri    = u[2] - u[1]
  rj    = u[3] - u[1]
  a     = 50.0 * (3pi/4 - pi/2) / (sqrt(1 - cos(3pi/4)^2) * norm(ri) * norm(rj))
  Z[2] .= a .* (rj .- (ri .* (dot(ri, rj) / dot(ri,ri))))
  Z[3] .= a .* (ri .- (rj .* (dot(ri, rj) / dot(rj,rj))))
  Z[1] .= -1 .* (Z[2] .+ Z[3])

  E1    = JMD._harmonicBondAngle!(F, u, 2, 1, 3, 50.0, pi/2)
  E2, f = JMD._harmonicBondAngle(u[2]-u[1], u[3]-u[1], 50.0, pi/2)
  @test isapprox(E1, E; atol=1e-8)
  @test isapprox(E2, E; atol=1e-8)
  @test isapprox(F[1], Z[1]; atol=1e-8)
  @test isapprox(F[2], Z[2]; atol=1e-8)
  @test isapprox(F[3], Z[3]; atol=1e-8)


end


