using Aqua
using Test
using JMD
using LinearAlgebra

@testset "Aqua.jl" begin
  Aqua.test_all(JMD)
end

@testset "IO" begin
  include("./io_test.jl")
end

@testset "Helpers" begin
  include("./helpers_test.jl")
end

@testset "Cells" begin
  include("./md/cells_test.jl")
end

@testset "Potentials" begin

  @testset "HGNN" begin
    include("./potentials/hgnn_test.jl")
  end

  @testset "MvHff" begin
    include("./potentials/cocoff_test.jl")
  end

end

@testset "Math Toolkit" begin

  @testset "Peak Finding" begin
    include("./mathtk/peakFinding_test.jl")
  end

end

@testset "Structural" begin

  @testset "Distributions" begin
    include("./structural/distributions_test.jl")
  end

end