__precompile__()
"""
Julia Molecular Dynamics

A very minimal package for me to do MD.

author: Brian C. Ferrari
"""

module JMD

using StaticArrays

include("./io.jl")
include("./md/bodies.jl")
include("./md/potentials/COCOff.jl")
include("./md/simulation.jl")

end # module
