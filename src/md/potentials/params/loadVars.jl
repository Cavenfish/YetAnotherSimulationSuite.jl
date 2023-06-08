f = jldopen(joinpath(@__DIR__, "hgnn.jld2"))
const hgnnPairVars = f["pairVars"]
const hgnnMolVars  = f["molVars"]
close(f)
