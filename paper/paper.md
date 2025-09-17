---
title: "YetAnotherSimulationSuite.jl: An Atomic Simulation Suite in Julia"
tags:
    - Julia
    - Atomic Simulations
    - Molecular Dynamics

authors:
    - name: "Brian C. Ferrari"
      orcid: 0000-0002-7416-8629
      affiliation: 1

affiliations:
    - name: "Leiden Institute of Chemistry, Leiden University, Leiden 2300 RA, The Netherlands"
      index: 1

date: 
bibliography: paper.bib
---

# Summary
`YetAnotherSimulationSuite.jl` (YASS) is a simulation suite
for atomic simulations in Julia [@julia]. YASS aims to be
highly performant, memory efficient, flexible, and extensible.
YASS was developed with a focus on making it incredibly easy
to add or customize anything within the package, making niche
research methods more accessible to the average user. 

YASS features a wide variety of capabilities useful for
molecular simulations:

- Reading and writing over 20 different file types
- Basic manipulation of atomic structures
- Geometry and cell optimizations
- Harmonic frequency calculations
- Velocity autocorrelation
- Molecular dynamics simulations
- Radial and angular distribution functions

# Statement of Need

Within the field of atomic simulations, there exists a vast
number of software packages for performing these simulations.
There are those that focus on ease of use and flexibility (i.e.,
ASE [@ase]), and those that focus on speed and memory efficiency
(i.e., LAMMPS[@lammps], GROMACS[@gromacs1; @gromacs2; @gromacs3;
@gromacs4], JaxMD[@jaxmd]). In an effort to offer both ease of use
and speed, ASE offers interfaces to many of the performance-
focused simulation suites. However, using these interfaces limits
the flexibility originally offered by ASE and can reduce the
memory efficiency offered by the faster simulation suite. The
Julia programming language offers the perfect solution to this
problem, as it can be as performant as C++ with the simplicity
of Python.
 
In theoretical chemistry, it is quite common for two packages with
similar functionality to flourish; for instance, both LAMMPS and
GROMACS offer similar functionality, but users tend to prefer the
interface of one over the other. This helps ensure that for all
types of users there is an interface that feels more natural for
them. Currently, `Molly.jl` [@molly] and `NQCDynamics.jl`
[@nqcdynamics] are the only packages in the Julia ecosystem that
perform molecular dynamics (MD) simulations. The latter focuses
on the more niche topic of nonadiabatic quantum classical dynamics
(NQCD), leaving only the former as an option for users wanting a
general atomic simulation suite. This creates a problem for users
who dislike the `Molly.jl` interface but want to perform MD in
Julia. YASS solves this problem by offering users an alternative
interface for atomic simulations.

# Examples

Two simple examples are shown here to illustrate how YASS can
be used for vibrational frequency analysis of molecular systems.
The examples focus only on this topic, as the process also
utilizes many additional YASS features. The first example below
shows the steps required to calculate the harmonic frequencies of
a water molecule.

```julia
using Optim
using YetAnotherSimulationSuite

# Read initial structure
molecule = readSystem("h2o.xyz")

# Run geometry optimization
optimized = opt(TIP4Pf(), LBFGS(), molecule)

# Save optimized structure
write("optimized.xyz", optimized)

# Calculate frequencies and normal modes
freqs, modes = getHarmonicFreqs(TIP4Pf(), optimized)
```

The second example, shown below, calculates the vibrational
density of states (VDOS) using the velocity autocorrelation
function (VACF). In this example, the `water_at_250K.xyz` file
should contain a cell with water at 250 Kelvin. YASS can also be
used to generate this initial configuration, but for simplicity,
it is not shown here. 

```julia
using YetAnotherSimulationSuite

# Read initial structure
water = readSystem("water_at_250K.xyz")

# Create NVE ensemble
ensemble = NVE(water)

# Run 10 picosecond simulation with 0.1 fs timestep
traj = run(TIP4Pf(), water, 10u"ps", 0.1u"fs", ensemble)

# Extract velocities and masses
vel, mas = getVelMas(traj)

# Configure VACF calculation
inp = vacfInps(
    vel,       # Velocity trajectories
    mas,       # Atomic masses
    1e16u"Hz", # Sampling frequency (1/fs = 1e15 Hz)
    true,      # Normalize VACF
    Hann,      # Window function
    4,         # FFT padding factor
    true       # Mirror the data
)

# Calculate VDOS
out = VDOS(inp)
```

More examples and detailed explanations can be found in the YASS
[documentation](https://cavenfish.github.io/YetAnotherSimulationSuite.jl/dev/).

# Dependencies

YASS relies on several packages and libraries that deserve to be
credited.

- `Chemfiles` is used for I/O operations.
- `Optim.jl`[@optim] is used for geometry and cell optimization.
- `OrdinaryDiffEq.jl` [@diffeq] is used as the integrator for MD simulations.
- `FiniteDifferences.jl` is used for computing Jacobians for harmonic frequency analysis.
- `FFTW` [@fftw] is used for fast Fourier transforms.
- `StaticArrays.jl` is used for memory efficiency.
- `JLD2.jl` is used to store complex data structures.
- `Clustering.jl` is used for molecular identification.

# Acknowledgements

BCF thanks Katie Slavicinska for creating the YASS logo. 

# References