__precompile__()
"""
Julia Molecular Dynamics

A very minimal package for me to do MD.

author: Brian C. Ferrari
"""

module JMD

using StaticArrays
using OrdinaryDiffEq

include("./io.jl")
include("./md/bodies.jl")
include("./md/potentials/COCOff.jl")
include("./md/potentials/HGNN.jl")
include("./md/simulation.jl")

end # module
