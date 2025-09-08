# Bodies

There are two types of simulation objects within YASS, the Particle and Cell objects. A Cell object is self-contained and is passed as-is to simulation functions, whereas the Particle object is only a single Particle. For simulations with more than one Particle, you must use a vector of Particle objects.

### Particle

In YASS you have an Particle object, which contains the position, velocity, mass and symbol of the atom. As mentioned above, you will typically be working with a vector of atoms. You can create this vector of atoms manually or read them in from an xyz file.

```julia
using YASS

# Read it in from an xyz file
bdys = readSystem("./myFile.xyz")

# Or make it yourself
bdys::Vector{YASS.MyAtoms} = [
  Particle([ 0.00,  0.22, 0.0], zeros(3), 15.999, "O"),
  Particle([ 0.75, -0.36, 0.0], zeros(3),  1.000, "H"),
  Particle([-0.75, -0.36, 0.0], zeros(3),  1.000, "H")
]

# Write it to an xyz file
write("./h2o.xyz", bdys)
```

Within YASS the Particle struct is mutable, this is to allow on-the-fly swapping of atomic masses. 

```julia
using YASS

# Read it in from an xyz file
bdys = readSystem("./myFile.xyz")

# Change the mass of the first atom
bdys[1].m = 2.00
```

YASS also contains some auxiliary functions that can apply changes to Particle objects.

```julia
using YASS

# Read it in from an xyz file
bdys = readSystem("./myFile.xyz")

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

The Cell object within YASS holds the lattice, scaled positions, velocities, masses, symbols, and some PBC criteria. Making the Cell object completely from sctratch is possible but not typical. 

```julia
using YASS

# Read it in from an xyz file
cell = readSystem("./myFile.xyz")

# Or make it yourself
bdys::Vector{YASS.MyAtoms} = [
  Particle([ 0.00,  0.22, 0.0], zeros(3), 15.999, 'O'),
  Particle([ 0.75, -0.36, 0.0], zeros(3),  1.000, 'H'),
  Particle([-0.75, -0.36, 0.0], zeros(3),  1.000, 'H')
]
lat  = [5 0 0; 0 5 0; 0 0 5]
cell = makeCell(bdys, lat)

# Write it to an xyz file
write("./myCell.xyz", cell)
```

YASS also contians some basic auxiliary functions for cells. 

```julia
using YASS

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