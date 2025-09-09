# Molecular Dynamics Simulations

`YASS.jl` provides functionality for classical molecular dynamics simulations in different ensembles. This guide explains how to set up and run MD simulations.

## Basic Usage

The simplest way to run an MD simulation is in the NVE (microcanonical) ensemble:

```julia
using YASS

# Read initial structure
water = readSystem("water.xyz")

# Create NVE ensemble
ensemble = NVE()

# Run 5 picosecond simulation with 0.1 fs timestep
traj = run(TIP4Pf(), water, (0.0, 5ps), 0.1fs, ensemble)
```

The `run` function takes the following arguments:
- Calculator (force field)
- Initial structure
- Time span tuple (start, end)
- Time step
- Ensemble

## Available Ensembles

### NVE Ensemble
The NVE ensemble maintains constant number of particles (N), volume (V), and energy (E):

```julia
ensemble = NVE()
```

### NVT Ensemble 
The NVT ensemble maintains constant temperature using a thermostat. `YASS.jl` supports several thermostats:

```julia
# Berendsen thermostat at 300K with 50fs coupling time
ensemble = Berendsen(300.0, 50fs) |> NVT

# Nosé-Hoover thermostat at 300K
ensemble = NoseHoover(300.0) |> NVT 

# Langevin thermostat at 300K with γ=1.0
ensemble = Langevin(300.0, 1.0) |> NVT
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
traj = run(TIP4Pf(), water, (0.0, 1ns), 1.0fs, ensemble; split=10)
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

## Periodic Boundary Conditions

For periodic systems, read the structure as a cell:

```julia
# Read periodic cell
cell = readSystem("crystal.xyz")

# Run simulation with PBC
traj = run(calc, cell, (0.0, 10ps), 1.0fs, NVE(cell))
```