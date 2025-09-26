@testset "Hit and Stick Call Test" begin
  calc = TIP4Pf()
  core = joinpath(@__DIR__, "../testingFiles/xyzFiles/h2o-dimer.xyz")
  mol  = joinpath(@__DIR__, "../testingFiles/xyzFiles/h2o.xyz")
  out  = joinpath(@__DIR__, "../testingFiles/xyzFiles/out.xyz")

  inp = HnS(
    30,
    5u"ps",
    5u"ps",
    8.6e-4u"eV",
    core,
    mol,
    out,
    1u"fs",
    CVR(25.0u"K", 100u"fs", calc)
  )

  hitAndStick(calc, inp)

  bdys = readSystem(out)
  mols = getMols(bdys, 1.2)

  @test length(mols) == 10
  
  rm(out)

end
