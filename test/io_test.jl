
function bdysTesting(file, N, s, m, r, v)
  bdys = readXyz(file)

  @test length(bdys) == N
  @test [i.s for i in bdys] == s
  @test [i.m for i in bdys] == m
  @test [i.r for i in bdys] == r
  @test [i.v for i in bdys] == v
end

function cellTesting(file, N, lat)
  cell = readCell(file)

  @test length(cell.masses) == N
  @test cell.lattice == lat
end

@testset "Reading files" begin

  bdysTesting(
    "./testingFiles/xyzFiles/co.xyz",
    2, ['C', 'O'], [12.011, 15.999], 
    [[0.0,0.0,0.0],[1.1,0.0,0.0]], 
    [[0.0,0.0,0.0],[0.0,0.0,0.0]]
  )

  bdysTesting(
    "./testingFiles/xyzFiles/co2.xyz",
    3, ['C', 'O', 'O'], [12.011, 15.999, 15.999], 
    [[0.0,0.0,0.0],[-1.12,0.0,0.0],[1.12,0.0,0.0]], 
    [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
  )

  bdysTesting(
    "./testingFiles/xyzFiles/ch4.xyz",
    5, ['C', 'H', 'H', 'H', 'H'], [12.011, 1.007, 1.007, 1.007, 1.007], 
    [[0.0,0.0,0.0],[0.0,0.0,1.0],[1.0,0.0,0.0],[-0.5,-0.8,-0.3],[-0.5,0.8,-0.3]], 
    [[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]]
  )

  cellTesting(
    "./testingFiles/xyzFiles/iceIh.xyz",
    2304, [30.9 0.0 0.0; 0.0 26.8 0.0; 0.0 0.0 29.1]
  )
  
end