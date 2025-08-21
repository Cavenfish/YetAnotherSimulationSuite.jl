
function bdysTesting(f, N, s, m, r, v)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)

  @test length(bdys) == N
  @test [i.s for i in bdys] == s
  @test [i.m for i in bdys] == m
  @test [i.r for i in bdys] == r
  @test [i.v for i in bdys] == v

  write("./tmp.xyz", bdys)
  new = readSystem("./tmp.xyz")
  rm("./tmp.xyz")

  @test new[end].r == bdys[end].r
  @test new[end].v == bdys[end].v
  @test new[end].m == bdys[end].m
  @test new[end].s == bdys[end].s
end

function cellTesting(f, N, lat)
  file = joinpath(@__DIR__, f)
  cell = readSystem(file)

  @test length(cell.masses) == N
  @test cell.lattice == lat

  write("./tmp.xyz", cell)
  new = readSystem("./tmp.xyz")
  rm("./tmp.xyz")

  @test cell.lattice == new.lattice
  @test cell.scaled_pos[123] ≈ new.scaled_pos[123] atol = 1e-5
end

function trajTesting(f, N)
  file = joinpath(@__DIR__, f)
  traj = readSystem(file)

  @test length(traj.images) == N

  write("./tmp.xyz", traj)
  new = readSystem("./tmp.xyz")
  rm("./tmp.xyz")

  @test length(new.images) == N
end

@testset "YASS Style File IO" begin

  bdysTesting(
    "testingFiles/xyzFiles/co.xyz",
    2, ["C", "O"], [12.011, 15.999], 
    [[0.0,0.0,0.0],[1.1,0.0,0.0]], 
    [[0.0,0.0,0.0],[0.0,0.0,0.0]]
  )

  bdysTesting(
    "testingFiles/xyzFiles/co2.xyz",
    3, ["C", "O", "O"], [12.011, 15.999, 15.999], 
    [[0.0,0.0,0.0],[-1.12,0.0,0.0],[1.12,0.0,0.0]], 
    [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
  )

  bdysTesting(
    "testingFiles/xyzFiles/ch4.xyz",
    5, ["C", "H", "H", "H", "H"], [12.011, 1.008, 1.008, 1.008, 1.008], 
    [[0.0,0.0,0.0],[0.0,0.0,1.0],[1.0,0.0,0.0],[-0.5,-0.8,-0.3],[-0.5,0.8,-0.3]], 
    [[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]]
  )

  cellTesting(
    "testingFiles/xyzFiles/iceIh.xyz",
    2304, [30.9 0.0 0.0; 0.0 26.8 0.0; 0.0 0.0 29.1]
  )

  trajTesting(
    "testingFiles/mdData/traj.xyz",
    51
  )
  
end

@testset "Read ASE xyz" begin
  file = joinpath(@__DIR__, "testingFiles/xyzFiles/co_ase.xyz")
  bdys = readSystem(file)

  @test bdys[1].s == "C"
  @test bdys[1].m == 12.011
  @test bdys[2].s == "O"
  @test bdys[2].m == 15.999
  @test bdys[2].v == [0.0, 0.0, 0.0]
  @test isapprox(bdys[1].r, [-0.07311636, -0.11425678, -0.02537179], atol=1e-10)

  file = joinpath(@__DIR__, "testingFiles/xyzFiles/Ih_ase.xyz")
  cell = readSystem(file)
  lat  = [7.82 0.0 0.0; -3.91 6.77 0.0; 0.0 0.0 7.36]

  @test isapprox(cell.lattice, lat, atol=1e-2)
end