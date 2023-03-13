__precompile__()
"""
Julia Molecular Dynamics

A very minimal package for me to do MD.

author: Brian C. Ferrari
"""

module JMD

using NBodySimulator
using StaticArrays

include("./io.jl")
include("./bodies.jl")
include("./COCOff.jl")

export readASExyz, COCOdyn

end # module
