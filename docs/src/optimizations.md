# Optimizations

JMD uses [Optim.jl]() for optimizations. You can optimize the geometry of a system or the lattice vectors of a cell. 

```julia
using JMD

# Read in bdys
bdys = readXyz("./myFile.xyz")

# Run a geometry optimization
new  = opt(TIP4P, JMD.LBFGS(), bdys)
```