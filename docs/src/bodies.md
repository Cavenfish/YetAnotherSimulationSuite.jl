# Bodies

There are two types of simulation objects within JMD, the Atom and Cell objects. A Cell object is self-contained and is passed as-is to simulation functions, whereas the Atom object is only a single Atom. For simulations with more than one Atom, you must use a vector of Atom objects.

### Atom

```julia
using JMD

# Read it in from an xyz file
bdys = readXyz("./myFile.xyz")

# Or make it yourself
bdys = [
  Atom([ 0.00,  0.22, 0.0], zeros(3), 15.999, 'O'),
  Atom([ 0.75, -0.36, 0.0], zeros(3),  1.000, 'H'),
  Atom([-0.75, -0.36, 0.0], zeros(3),  1.000, 'H')
]

# Write it to an xyz file
writeXyz("./h2o.xyz", bdys)
```

### Cell