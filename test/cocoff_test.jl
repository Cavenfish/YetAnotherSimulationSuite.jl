using StaticArrays
using Test

include("../src/JMD.jl")

@testset "CO-CO FF Test" begin

  @testset "10 CO Molecules" begin
    bdys = JMD.readASExyz("10co.xyz")
    pos   = [SVector{3}(i.r) for i in bdys]
    mas   = [i.m for i in bdys]   
    pars, mols = JMD.getPairs(bdys)
    p   = JMD.NVEsimu(bdys, pars, mols, [], [], mas)
    JMD.COCO(zero(pos),zero(pos),pos,p,0)

    #Energy Tests
    @test p.energy[1]["Morse"] ≈  1.7002672151988918e-05 atol=1e-15
    @test p.energy[1]["Exch"]  ≈  0.285819645468938      atol=1e-12
    @test p.energy[1]["Disp"]  ≈ -0.5667142867222983     atol=1e-12
    @test p.energy[1]["Coul"]  ≈ -0.0786977801615315     atol=1e-12

    #Force Tests
    @test p.forces[1]["Morse"][6] ≈ [-0.00087197, -0.00224147,  0.00105403] atol=1e-8
    @test p.forces[1]["Exch"][ 9]  ≈ [-0.04894114, -0.06273592, -0.04394393] atol=1e-8
    @test p.forces[1]["Disp"][13]  ≈ [ 0.02416868, -0.04554598, -0.00693859] atol=1e-8
    @test p.forces[1]["Coul"][16]  ≈ [-0.00378278, -0.01104432, -0.0009977 ] atol=1e-8

  end

end
