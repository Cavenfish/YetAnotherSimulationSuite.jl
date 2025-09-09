# Simulation Bodies

`YASS.jl` provides two main types of simulation objects: `Particle` and `Cell`. Each serves different purposes in molecular simulations:

- `Particle`: Represents individual atoms or particles
- `Cell`: Represents periodic systems like crystals

## Particle Objects

The `Particle` type is the fundamental building block for molecular simulations. Each particle has:

- Position vector (`r`)
- Velocity vector (`v`) 
- Mass (`m`)
- Chemical symbol (`s`)

### Creating Particles

You can create particles in several ways:

```julia
using YASS

# Read from XYZ file
atoms = readSystem("molecule.xyz")

# Create manually
water::Vector{YASS.MyAtoms} = [
    Particle([0.000,  0.000, 0.000], zeros(3), 15.999, "O"),
    Particle([0.757,  0.586, 0.000], zeros(3),  1.008, "H"),
    Particle([0.757, -0.586, 0.000], zeros(3),  1.008, "H")
]

# Save to file
write("water.xyz", water)
```

### Modifying Particles

The `Particle` type is mutable, allowing modifications after creation:

```julia
# Modify position
atoms[1].r .= [1.0, 0.0, 0.0]

# Change velocity
atoms[1].v .= [0.1, 0.0, 0.0]

# Update mass (e.g., for isotope studies)
atoms[1].m = 18.015  # Change to heavy water

# Change chemical symbol
atoms[1].s = "D"     # Deuterium
```

### Particle Manipulation Functions

YASS provides several utility functions for working with particles:

```julia
# Translate all particles
translate!(atoms, [1.0, 0.0, 0.0])

# Center particles at origin
centerBdys!(atoms)

# Swap positions of two particles
swapAtoms!(atoms, 1, 2)

# Change isotopes for specific atoms
swapIso!(atoms, [1,2], [2.014, 2.014])  # Convert H to D

# Get center of mass
com = CoM(atoms)

# Get center of mass velocity
vcom = vCoM(atoms)

# Remove center of mass motion
zeroVCoM!(atoms)
```

## Cell Objects

The `Cell` type represents periodic systems and contains:

- Lattice matrix
- Scaled positions (fractional coordinates)
- Velocities
- Masses
- Chemical symbols  
- Periodic boundary conditions
- Neighbor counts

### Creating Cells

```julia
using YASS

# Read from file with lattice information
cell = readSystem("crystal.xyz")

# Create from atoms and lattice
atoms = [
    Particle([0.0, 0.0, 0.0], zeros(3), 22.990, "Na"),
    Particle([0.5, 0.5, 0.5], zeros(3), 35.450, "Cl")
]
lattice = [
    5.0 0.0 0.0
    0.0 5.0 0.0
    0.0 0.0 5.0
]
cell = makeCell(atoms, lattice)

# Specify periodic boundary conditions
cell = makeCell(atoms, lattice, 
    PBC=[true, true, true],  # Periodic in all directions
    NC=[1,1,1]               # Neighbor cells to consider
)

# Save cell to file
write("nacl.xyz", cell)
```

### Cell Operations

`YASS.jl` provides various functions for cell manipulation:

```julia
# Wrap atoms back into primary cell
wrap!(cell)

# Center atoms in cell
center!(cell)

# Create supercell
transform = [2 0 0; 0 2 0; 0 0 2]  # 2x2x2 supercell
super = makeSuperCell(cell, transform)

# Get primitive cell (with symmetry precision)
primitive = getPrimitiveCell(cell, 1e-5)

# Convert between cell and atoms
atoms = makeBdys(cell)        # Cell -> Atoms
cell = makeCell(atoms, lat)   # Atoms -> Cell

# Get Cartesian positions
positions = getPos(cell)

# Get cell volume
volume = getVolume(cell)

# Reorder atoms
order = sortperm([atom.m for atom in atoms])
reorder!(cell, order)
```