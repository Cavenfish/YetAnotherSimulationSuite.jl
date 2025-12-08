# Molecular Dynamics Simulations

YASS provides functionality for classical molecular dynamics simulations in different ensembles. This guide explains how to set up and run MD simulations.

## Basic Usage

The simplest way to run an MD simulation is in the NVE (microcanonical) ensemble:

```julia
using YetAnotherSimulationSuite

# Read initial structure
water = readSystem("water.xyz")

# Create NVE ensemble
ensemble = NVE()

# Run 5 picosecond simulation with 0.1 fs timestep
traj = run(TIP4Pf(), water, 5u"ps", 0.1u"fs", ensemble)

# You can also specify the start time
traj = run(TIP4Pf(), water, (5u"ps", 10u"ps"), 0.1u"fs", ensemble)
```

The `run` function takes the following arguments:
- Calculator (force field)
- Initial structure
- Time or time span tuple (start, end)
- Time step
- Ensemble

## Available Ensembles

### NVE Ensemble
The NVE ensemble maintains constant number of particles (N), volume (V), and energy (E):

```julia
ensemble = NVE()
```

### NVT Ensemble 
The NVT ensemble maintains constant temperature using a thermostat. YASS supports several thermostats:

```julia
# Thermostats need unit information from calc
calc = TIP4Pf()

# Berendsen thermostat at 300K with 50fs coupling time
ensemble = Berendsen(300.0u"K", 50u"fs", calc) |> NVT

# Cannonical velocity rescaling thermostat at 300K
ensemble = CVR(300.0u"K", 100u"fs", calc) |> NVT 

# Langevin thermostat at 300K with 10fs coupling time
ensemble = Langevin(300.0u"K", 10u"fs") |> NVT
```

## Working with Trajectories

The `run` function returns a `Traj` object containing the simulation data:

```julia
# Access trajectory information
println("Number of frames: ", length(traj))
println("Atomic symbols: ", traj.symbols)
println("Atomic masses: ", traj.masses)

# Get positions from first frame
pos = traj.images[1].pos

# Get velocities from last frame 
vel = traj.images[end].vel

# Get temperature and energy arrays
temps = [img.temp for img in traj.images]
energies = [img.energy for img in traj.images]
```

## Long Simulations

For longer simulations, you can split them into segments to save memory:

```julia
# Run 1 nanosecond simulation split into 10 segments
traj = run(TIP4Pf(), water, 1u"ns", 1.0u"fs", ensemble; split=10)
```

## Saving Trajectories

Trajectories can be saved to disk in various formats:

```julia
# Save as XYZ file
write("trajectory.xyz", traj)

# Save as JLD file (binary format)
using JLD2
@save "trajectory.jld2" traj
```

## Periodic Boundary Conditions (PBC)

For periodic systems, read the structure as a cell:

```julia
# Read periodic cell
cell = readSystem("crystal.xyz")

# Run simulation with PBC
traj = run(calc, cell, 10u"ps", 1.0u"fs", NVE(cell))
```