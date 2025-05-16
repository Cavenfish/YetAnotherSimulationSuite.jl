
function staticTest(EoM, f)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)

  getPotEnergy(EoM, bdys)
  getForces(EoM, bdys)

  true
end

function dynamicTest(EoM, f)
  file = joinpath(@__DIR__, f)
  bdys = readSystem(file)

  runMD(EoM, bdys, (0.0, 10fs), 1fs)

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
  @test staticTest(MvHff, co)
  @test dynamicTest(MvHff, co)

  # HGNN
  @test staticTest(HGNN, co)
  @test dynamicTest(HGNN, co)

  # TIP4P
  @test staticTest(TIP4P, h2o)
  @test dynamicTest(TIP4P, h2o)

  # SPCF
  @test staticTest(SPCF, h2o)
  @test dynamicTest(SPCF, h2o)

end