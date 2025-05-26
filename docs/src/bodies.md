# Bodies

There are two types of simulation objects within JMD, the Atom and Cell objects. A Cell object is self-contained and is passed as-is to simulation functions, whereas the Atom object is only a single Atom. For simulations with more than one Atom, you must use a vector of Atom objects.

### Atom

In JMD you have an Atom object, which contains the position, velocity, mass and symbol of the atom. As mentioned above, you will typically be working with a vector of Atoms. You can create this vector of Atoms manually or read them in from an xyz file.

```julia
using JMD

# Read it in from an xyz file
bdys = readSystem("./myFile.xyz")

# Or make it yourself
bdys = [
  Particle([ 0.00,  0.22, 0.0], zeros(3), 15.999, 'O'),
  Particle([ 0.75, -0.36, 0.0], zeros(3),  1.000, 'H'),
  Particle([-0.75, -0.36, 0.0], zeros(3),  1.000, 'H')
]

# Write it to an xyz file
write("./h2o.xyz", bdys)
```

Within JMD the Atom struct is mutable, this is to allow on-the-fly swapping of atomic masses. 

```julia
using JMD

# Read it in from an xyz file
bdys = readSystem("./myFile.xyz")

# Change the mass of the first Atom
bdys[1].m = 2.00
```

JMD also contains some auxiliary functions that can apply changes to Atom objects.

```julia
using JMD

# Read it in from an xyz file
bdys = readXyz("./myFile.xyz")

# Uniformally translate all atoms
v = [1.0, 0.0, 0.0]
translateBdys!(bdys, v)

# Center the atoms about the center-of-mass
centerBdys!(bdys)

# Swap the positions of atoms 1 and 2
swapAtoms!(bdys, 1, 2)

# Swap the mass of atoms 1, and 2 with 4.00 and 5.00
swapIso!(bdys, [1,2], [4.0, 5.0])
```

### Cell

The Cell object within JMD holds the lattice, scaled positions, velocities, masses, symbols, and some PBC criteria. Making the Cell object completely from sctratch is possible but not typical. 

```julia
using JMD

# Read it in from an xyz file
cell = readSystem("./myFile.xyz")

# Or make it yourself
bdys = [
  Particle([ 0.00,  0.22, 0.0], zeros(3), 15.999, 'O'),
  Particle([ 0.75, -0.36, 0.0], zeros(3),  1.000, 'H'),
  Particle([-0.75, -0.36, 0.0], zeros(3),  1.000, 'H')
]
lat  = [5 0 0; 0 5 0; 0 0 5]
cell = makeCell(bdys, lat)

# Write it to an xyz file
write("./myCell.xyz", cell)
```

JMD also contians some basic auxiliary functions for cells. 

```julia
using JMD

# Read it in from an xyz file
cell = readSystem("./myFile.xyz")

# Wrap atoms outside cell back into cell
wrap!(cell)

# Center atoms at cell center
center!(cell)

# Make a supercell
T     = [2 0 0; 0 2 0; 0 0 2]
super = makeSuperCell(cell, T)

# Get primitive cell
prim  = getPrimitiveCell(cell, 1e-5)# 1e-5 is symprec

# Get vector of Atoms from cell
bdys  = getBdys(cell)
```