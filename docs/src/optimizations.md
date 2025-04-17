# Optimizations

JMD uses [Optim.jl]() for optimizations. You can optimize the geometry of a system or the lattice vectors of a cell. 

```julia
using JMD

# Read in bdys
bdys = readXyz("./myFile.xyz")

# Run a geometry optimization
new  = opt(TIP4P, JMD.LBFGS(), bdys)
```

The optimization algorithms from Optim.jl are not re-exported by JMD, for this reason you must either use `JMD.LBFGS()` (as shown above) or load Optim.jl within your script (as shown below).

```julia
# Alternatively
using JMD
using Optim

# Read in bdys
bdys = readXyz("./myFile.xyz")

# Run a geometry optimization
new  = opt(TIP4P, LBFGS(), bdys)
```

