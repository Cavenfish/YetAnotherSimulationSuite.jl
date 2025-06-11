
function staticTest(calc, f)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)

  getPotEnergy(calc, bdys)
  getForces(calc, bdys)

  true
end

function dynamicTest(calc, f)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)

  run(calc, bdys, (0.0, 10fs), 1fs, NVE())

  true
end

function thermoTest(calc, thermo, f)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)

  run(calc, bdys, (0.0, 10fs), 1fs, NVT(thermo); split=2)

  true
end

# This collection of tests simply calls a static 
# and dynamics potential call. It does not test 
# if the potentials perform correctly, rather it 
# checks if any code changes breaks their calls. 
@testset "Test Potential Function Calls" begin
  co  = "../testingFiles/xyzFiles/co-dimer.xyz"
  h2o = "../testingFiles/xyzFiles/h2o-dimer.xyz"
  
  # MvHff
  @test staticTest(MvHff(), co)
  @test dynamicTest(MvHff(), co)

  # HGNN
  @test staticTest(HGNN(), co)
  @test dynamicTest(HGNN(), co)

  # TIP4P
  @test staticTest(TIP4Pf(), h2o)
  @test dynamicTest(TIP4Pf(), h2o)

  # SPCF
  @test staticTest(SPCF(), h2o)
  @test dynamicTest(SPCF(), h2o)

end

@testset "Test Thermostat Function Calls" begin
  h2o = "../testingFiles/xyzFiles/h2o-dimer.xyz"

  # Berendsen
  @test thermoTest(TIP4Pf(), Berendsen(75.0, 100fs), h2o)

  # Langevin
  @test thermoTest(TIP4Pf(), Langevin(75.0, 100fs), h2o)

  # CVR
  @test thermoTest(TIP4Pf(), CVR(75.0, 100fs), h2o)
end