# Vibrational Analysis

`YASS.jl` provides multiple methods for analyzing vibrational properties of molecular systems:
1. Harmonic frequency analysis
2. Velocity autocorrelation function (VACF)

## Harmonic Frequencies

The harmonic approximation calculates vibrational frequencies by diagonalizing the mass-weighted Hessian matrix:

```julia
using YASS

# Read molecular structure
molecule = readSystem("water.xyz")

# Calculate frequencies and normal modes
freqs, modes = getHarmonicFreqs(TIP4Pf(), molecule)
```

The outputs are:
- `freqs`: Vector of vibrational frequencies (in cm$^{-1}$)
- `modes`: Matrix where each column is a normal mode eigenvector

### Analyzing Normal Modes

You can visualize and analyze individual modes:

```julia
# Get the first normal mode
mode1 = modes[:,1]

# Animate a specific mode
animateMode(molecule, mode1, "mode1.xyz", c=0.5)  # c controls amplitude

# Calculate potential energy surface along mode
pes = getModePES(TIP4Pf(), molecule, mode1)
```

### Mode Selection

For larger molecules, you can filter and analyze specific modes:

```julia
# Find modes in a frequency range
range = 3000:4000  # OH stretch region
idx = findall(f -> f in range, real.(freqs))
stretch_modes = modes[:,idx]

# Calculate mode participation ratios
pr = getPR(modes)  # Shows which atoms participate in each mode

# Get inverse participation ratio
ipr = getIPR(modes)
```

## Velocity Autocorrelation Function (VACF)

The VACF analyzes vibrational properties from MD trajectories:

```julia
using YASS

# Run MD simulation
molecule = readSystem("water.xyz")
traj = run(TIP4Pf(), molecule, (0.0, 10.0ps), 1.0fs, NVE())

# Extract velocities and masses
vel, mas = getVelMas(traj)

# Configure VACF calculation
inp = vacfInps(
    vel,           # Velocity trajectories
    mas,           # Atomic masses
    1e15,        # Sampling frequency (1/fs = 1e15 Hz)
    true,          # Normalize VACF
    Hann,          # Window function
    4,             # FFT padding factor
    true           # Mirror the data
)

# Calculate VDOS
out = VDOS(inp)
```

### VACF Components

The `vacfOut` structure contains:
- `out.c`: Raw velocity autocorrelation function
- `out.C`: Windowed/processed VACF
- `out.v`: Frequency axis (in cm$^{-1}$)
- `out.I`: Vibrational density of states

### Customizing the Analysis

Several parameters can be adjusted:

```julia
# Different window functions
inp = vacfInps(vel, mas, 1.0/fs, true, Welch, 4, true)   # Welch window
inp = vacfInps(vel, mas, 1.0/fs, true, HannM, 4, true)   # Modified Hann

# Increased padding for better frequency resolution
inp = vacfInps(vel, mas, 1.0/fs, true, Hann, 8, true)

# Without mirroring
inp = vacfInps(vel, mas, 1.0/fs, true, Hann, 4, false)
```

### Atom-Specific Analysis

You can analyze specific atoms or groups:

```julia
# Analyze only oxygen atoms
O_idx = findall(x -> x == "O", traj.symbols)
out_O = VDOS(inp, atms=O_idx)

# Analyze only hydrogen atoms
H_idx = findall(x -> x == "H", traj.symbols)
out_H = VDOS(inp, atms=H_idx)

# Compare spectra
using Plots
plot(out_O.v, out_O.I, label="Oxygen", alpha=0.6)
plot!(out_H.v, out_H.I, label="Hydrogen", alpha=0.6)
xlabel!("Wavenumber (cm-1)")
ylabel!("VDOS")
```

## Visualization

YASS provides several ways to visualize vibrational properties:

```julia
using Plots

# Plot VDOS spectrum
plot(out.v, out.I,
     xlabel="Wavenumber (cm⁻¹)",
     ylabel="VDOS",
     label="Total",
     linewidth=2)

# Plot raw VACF
plot(out.c,
     xlabel="Time",
     ylabel="VACF",
     label="Raw")

# Plot windowed VACF
plot(out.C,
     xlabel="Time",
     ylabel="VACF",
     label="Processed")
```

## Tips for Quality Results

1. For harmonic analysis:
   - Ensure structures are well-optimized
   - Use tight convergence criteria

2. For VACF analysis:
   - Use long enough trajectories (>10 ps)
   - Use an appropriate timestep (depends on mode frequency)
   - Ensure good energy conservation
   - Test different window functions
   - Adjust padding for desired resolution