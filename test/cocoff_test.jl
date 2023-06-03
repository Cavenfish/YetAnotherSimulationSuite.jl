using StaticArrays
using Test

include("../src/JMD.jl")

@testset "CO-CO FF Test" begin

  @testset "Single Molecule" begin
    c   = @SVector [-0.95671174, 2.73940642, 4.27069074]
    o   = @SVector [-0.2178486,  2.83265895, 3.42314641]
    co  = [c,o]
    res = JMD.COCOdyn(0,0,co,0,0)

    #Energy Tests
    @test res[1][1] ≈ 1.143575421025389e-07 atol=1e-15
    @test res[2][1] == 0.0
    @test res[3][1] == 0.0
    @test res[4][1] == 0.0

    #Force Tests
    @test res[1][2] ≈ [[ 0.00345525,  0.00043609, -0.00396349],
                       [-0.00345525, -0.00043609,  0.00396349]] atol=1e-8

    @test res[2][2] == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    @test res[3][2] == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    @test res[4][2] == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

  end #Single Molecule Testset

  @testset "10 CO Molecules" begin
    bdys = JMD.readASExyz("10co.xyz")
    pos  = [i.r for i in bdys]
    res  = JMD.COCOdyn(0,0,pos,0,0)

    #Energy Tests
    @test res[1][1] ≈  1.7002672151988918e-05 atol=1e-15
    @test res[2][1] ≈  0.285819645468938      atol=1e-12
    @test res[3][1] ≈ -0.5667142867222983     atol=1e-12
    @test res[4][1] ≈ -0.0786977801615315     atol=1e-12

    #Force Tests
    @test res[1][2][6 ] ≈ [-0.00087197, -0.00224147,  0.00105403] atol=1e-8
    @test res[2][2][9 ] ≈ [-0.04894114, -0.06273592, -0.04394393] atol=1e-8
    @test res[3][2][13] ≈ [ 0.02416868, -0.04554598, -0.00693859] atol=1e-8
    @test res[4][2][16] ≈ [-0.00378278, -0.01104432, -0.0009977 ] atol=1e-8

  end

end
