using Test, JMD, LinearAlgebra

@testset "IO" begin
  include("./io_test.jl")
end

@testset "Potentials" begin

  @testset "HGNN" begin
    include("./potentials/hgnn_test.jl")
  end

  @testset "MvHff" begin
    include("./potentials/cocoff_test.jl")
  end

end