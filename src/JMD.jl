__precompile__()
"""
Julia Molecular Dynamics

A very minimal package for me to do MD.

author: Brian C. Ferrari
"""

module JMD

using JLD2
using FFTW
using Optim
using StaticArrays
using Distributions
using OrdinaryDiffEq

const fs   = 0.09822694788464063 # 1fs in ASE time

include("./io.jl")
include("./helpers.jl")

include("./md/bodies.jl")
include("./md/potentials/COCOff.jl")
include("./md/potentials/HGNN.jl")
include("./md/simulation.jl")
include("./md/thermostats.jl")

include("./analysis/vacf.jl")
include("./analysis/vibrations.jl")
include("./analysis/optimizations.jl")


include("./md/potentials/params/loadVars.jl")
end # module
