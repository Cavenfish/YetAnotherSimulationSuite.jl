
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

  if isa(bdys, MyCell)
    run(calc, bdys, (0.0u"fs", 5u"fs"), 1u"fs", NVE(bdys))
    run(calc, bdys, 5u"fs", 1u"fs", NVE(bdys))
  else
    run(calc, bdys, (0.0u"fs", 5u"fs"), 1u"fs", NVE())
    run(calc, bdys, 5u"fs", 1u"fs", NVE())
  end

  true
end

function thermoTest(calc, thermo, f)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)

  run(calc, bdys, (0.0u"fs", 10u"fs"), 1u"fs", NVT(thermo); split=2)

  true
end

# This collection of tests simply calls a static 
# and dynamics potential call. It does not test 
# if the potentials perform correctly, rather it 
# checks if any code changes breaks their calls. 
@testset "Test Potential Function Calls" begin
  co  = "../testingFiles/xyzFiles/co-dimer.xyz"
  h2o = "../testingFiles/xyzFiles/h2o-dimer.xyz"
  ice = "../testingFiles/xyzFiles/iceIh_small.xyz"
  
  # MvHff
  @test staticTest(MvHff(), co)
  @test dynamicTest(MvHff(), co)

  # HGNN
  @test staticTest(HGNN(), co)
  @test dynamicTest(HGNN(), co)

  # TIP4P
  @test staticTest(TIP4Pf(), h2o)
  @test dynamicTest(TIP4Pf(), h2o)
  @test staticTest(TIP4Pf(), ice)
  @test dynamicTest(TIP4Pf(), ice)

  # SPCF
  @test staticTest(SPCF(), h2o)
  @test dynamicTest(SPCF(), h2o)
  @test staticTest(SPCF(), ice)
  @test dynamicTest(SPCF(), ice)

end

@testset "Test Thermostat Function Calls" begin
  h2o = "../testingFiles/xyzFiles/h2o-dimer.xyz"
  calc = TIP4Pf()

  # Berendsen
  @test thermoTest(calc, Berendsen(75.0u"K", 100u"fs", calc), h2o)

  # Langevin
  @test thermoTest(calc, Langevin(75.0u"K", 100u"fs", calc), h2o)

  # CVR
  @test thermoTest(calc, CVR(75.0u"K", 100u"fs", calc), h2o)
end