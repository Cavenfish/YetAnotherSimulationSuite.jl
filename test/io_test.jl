using SHA

function bdysTesting(f, N, s, m, r, v)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)
  og   = read(file) |> sha256

  @test length(bdys) == N
  @test [i.s for i in bdys] == s
  @test [i.m for i in bdys] == m
  @test [i.r for i in bdys] == r
  @test [i.v for i in bdys] == v

  write("./tmp.xyz", bdys)
  new = read("./tmp.xyz") |> sha256
  rm("./tmp.xyz")

  @test new == og
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

@testset "JMD Style File IO" begin

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
    5, ["C", "H", "H", "H", "H"], [12.011, 1.007, 1.007, 1.007, 1.007], 
    [[0.0,0.0,0.0],[0.0,0.0,1.0],[1.0,0.0,0.0],[-0.5,-0.8,-0.3],[-0.5,0.8,-0.3]], 
    [[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]]
  )

  cellTesting(
    "testingFiles/xyzFiles/iceIh.xyz",
    2304, [30.9 0.0 0.0; 0.0 26.8 0.0; 0.0 0.0 29.1]
  )
  
end

@testset "Read ASE xyz" begin
  file = joinpath(@__DIR__, "testingFiles/xyzFiles/co_ase.xyz")
  bdys = readSystem(file)

  @test bdys[1].s == "C"
  @test bdys[1].m == 12.011
  @test bdys[2].s == "O"
  @test bdys[2].m == 15.999
  @test bdys[1].r == [-0.07311636, -0.11425678, -0.02537179]
  @test bdys[2].v == [0.0, 0.0, 0.0]

  file = joinpath(@__DIR__, "testingFiles/xyzFiles/Ih_ase.xyz")
  
  bdys, cell = readSystem(file; getCell=true)

  @test cell == [
    7.82, 0.0, 0.0, 
    -3.9099999999999984, 6.772318657594311, 0.0,
    0.0, 0.0, 7.36
  ]
end