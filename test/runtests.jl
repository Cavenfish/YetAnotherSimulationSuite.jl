using Aqua
using Test
using Optim
using LinearAlgebra
using YetAnotherSimulationSuite

# Less typing
const YASS = YetAnotherSimulationSuite

@testset "Aqua.jl" begin
  Aqua.test_all(YetAnotherSimulationSuite)
end

@testset "IO" begin
  include("./io_test.jl")
end

@testset "Helpers" begin
  include("./helpers_test.jl")
end

@testset "Bodies" begin
  include("./md/bodies_test.jl")
end

@testset "Cells" begin
  include("./md/cells_test.jl")
end

@testset "Potentials" begin

  @testset "Intra Funcs" begin
    include("./potentials/funcs/intra_test.jl")
  end

  @testset "Cutoffs" begin
    include("./md/cutoff_test.jl")
  end

  @testset "Call Test" begin
    include("./potentials/calls_test.jl")
  end
  
  @testset "Test Specific Potentials" begin
    include("./potentials/cocoff_test.jl")
    include("./potentials/hgnn_test.jl")
    include("./potentials/lj_test.jl")
  end

end

@testset "Building" begin

  @testset "Hit and Stick" begin
    include("./building/hitAndStick_test.jl")
  end
  
end

@testset "Math Toolkit" begin

  @testset "Peak Finding" begin
    include("./mathtk/peakFinding_test.jl")
  end

  @testset "SavGol Filter" begin
    include("./mathtk/savitzkyGolay_test.jl")
  end

  @testset "Geometric Operations" begin
    include("./mathtk/geometric_operations_test.jl")
  end

  @testset "Stress Tensor" begin
    include("./mathtk/stress_test.jl")
  end

end

@testset "Structural" begin

  @testset "Distributions" begin
    include("./structural/distributions_test.jl")
  end

  @testset "Mols and Pairs" begin
    include("./structural/molsAndPairs_test.jl")
  end

end

@testset "Analysis" begin

  @testset "Energetics" begin
    include("./analysis/energetics_test.jl")
  end

  @testset "VACF" begin
    include("./analysis/vacf_test.jl")
  end

  @testset "Vibrations" begin
    include("./analysis/vibrations_test.jl")
  end

end