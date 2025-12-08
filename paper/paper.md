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

date: 22 September 2025
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

YASS offers users a similarly simple and easy-to-use interface
as ASE, while also offering significant speedups. JaxMD can also
offer a simple interface and high performance, but is limited
in the available optimization algorithms. YASS (through `Optim.jl`
[@optim]) offers a wide variety of optimization algorithms for
geometry and cell optimizations. Within the Julia ecosystem, 
`Molly.jl` [@molly] and `NQCDynamics.jl` [@nqcdynamics] are the 
dominant packages for performing molecular dynamics (MD) simulations.
However, YASS not only performs MD simulations but also
geometry/cell optimizations and harmonic frequency calculations.
A major shortcoming of the current state of YASS is a lack of 
support for parallelisation and GPU acceleration. Future
versions of YASS will aim to offer these features to
further improve performance.

# Examples

Two simple examples are shown here to illustrate how YASS can
be used for vibrational frequency analysis of molecular systems.
The examples focus only on this topic, as the process also
utilizes many additional YASS features. For both examples
the necessary xyz file can be found in the YASS 
[repo](https://github.com/Cavenfish/YetAnotherSimulationSuite.jl/tree/main/test/testingFiles/xyzFiles)
The first example below shows the steps required to calculate 
the harmonic frequencies of a water molecule.

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
function (VACF).  

```julia
using YetAnotherSimulationSuite

# Read initial structure
iceIh = readSystem("iceIh_small.xyz")

# Initialize calculator
calc = TIP4Pf()

# Create NVT ensemble with the
# Canonical Velocity Rescaling thermostat
# set to 50 Kelvin temp and 100 fs time constant
nvt_ensemble = NVT(iceIh, CVR(50.0u"K", 100u"fs", calc))

# Run 2 picosecond NVT simulation
nvt_traj = run(calc, iceIh, 2u"ps", 1u"fs", nvt_ensemble)

# Get last frame from NVT
bdys = makeBdys(nvt_traj, length(nvt_traj.images))

# Make it periodic
iceIh_50K = makeCell(bdys, iceIh.lattice)

# Create NVE ensemble
nve_ensemble = NVE(iceIh)

# Run 10 picosecond simulation with 0.1 fs timestep
nve_traj = run(calc, iceIh_50K, 10u"ps", 0.1u"fs", nve_ensemble)

# Extract velocities and masses
vel, mas = getVelMas(nve_traj)

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