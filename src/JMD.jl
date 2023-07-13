#__precompile__()
"""
Julia Molecular Dynamics

A very minimal package for me to do MD.

author: Brian C. Ferrari
"""

module JMD

using JLD2
using FFTW
using Optim
using MiniQhull
using DataFrames
using Statistics
using StaticArrays
using LinearAlgebra
using Distributions
using KernelDensity
using OrdinaryDiffEq
using FiniteDifferences

const fs   = 0.09822694788464063 # 1fs in ASE time
const ps   = 1000 * fs
const ns   = 1000 * ps
const kB   = 8.617333262e-5 # eV / K

include("./io.jl")
include("./helpers.jl")

include("./md/bodies.jl")
include("./md/potentials/COCOff.jl")
include("./md/potentials/HGNN.jl")
include("./md/simulation.jl")
include("./md/thermostats.jl")
include("./md/post-processing.jl")

include("./analysis/vacf.jl")
include("./analysis/vibrations.jl")
include("./analysis/optimizations.jl")
include("./analysis/structral.jl")

include("./mathtk/alphashape.jl")

include("./building/hitAndStick.jl")

include("./md/potentials/params/loadVars.jl")
end # module
