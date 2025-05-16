
# Read in some dummy bdys
co2 = readSystem("$(@__DIR__)/../testingFiles/xyzFiles/co2.xyz")
h2o = readSystem("$(@__DIR__)/../testingFiles/xyzFiles/h2o.xyz")

# Define a box lattice
lat = [10 0 0; 0 10 0; 0 0 10]

@testset "Test Basic Functions" begin
  cell = makeCell(co2, lat)
  spos = [inv(lat) * i.r for i in co2]
  pos  = [i.r for i in co2]
  posC = getPos(cell)

  @test isapprox(cell.scaled_pos, spos, atol=1e-10)
  @test isapprox(posC, pos, atol=1e-10)
  @test getVolume(cell) == 1000.0
  @test getMols(cell, 1.2) == [[1,2,3]]
  @test getPairs(cell)[1] |> isempty

end

@testset "Test Inplace Modifiers" begin
  cell = makeCell(h2o, lat)
  center!(cell)

  # check if molecule is centered
  tmp = cell.scaled_pos |> sum |> sum
  @test tmp / 9 == 0.5

  # save original atom spos
  tmp  = zeros(3)
  tmp .= cell.scaled_pos[1]

  # Move atom 1 a_vec away then wrap back
  cell.scaled_pos[1] .+= [10, 0, 0]
  wrap!(cell)

  @test isapprox(tmp, cell.scaled_pos[1], atol=1e-10)

  # delete hydrogens
  JMD.trim!(cell, 2:3)

  @test length(cell.scaled_pos) == 1
  @test length(cell.velocity)   == 1
  @test length(cell.masses)     == 1
  @test length(cell.symbols)    == 1

end

@testset "Test makeSuperCell" begin
  T     = [1 0 0; 0 2 0; 0 0 4]
  cell  = makeCell(h2o, lat)
  super = makeSuperCell(cell, T)

  @test super.lattice == [10.0 0 0; 0 20.0 0; 0 0 40.0]
  @test length(super.masses) == 3*2*4

  # grab some scaled coords
  z1 = super.scaled_pos[1][3]
  z2 = super.scaled_pos[4][3]
  z3 = super.scaled_pos[7][3]
  
  @test z2 - z1 == 1/4
  @test z3 - z1 == 2/4

  # save some super spos 
  check1  = zeros(3)
  check1 .= super.scaled_pos[1]
  check2  = zeros(3)
  check2 .= super.scaled_pos[8]
  check3  = zeros(3)
  check3 .= super.scaled_pos[12]

  # zero all values
  for i in super.scaled_pos
    i .= 0
  end

  # Hopefully recover all values
  makeSuperCell!(super, cell, T)

  @test check1 == super.scaled_pos[1]
  @test check2 == super.scaled_pos[8]
  @test check3 == super.scaled_pos[12]

end