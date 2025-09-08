# Vibrational Analysis

`YASS.jl` can calculate harmonic vibrational frequencies and eigenvectors. It is also possible to use a velocity autocorrelation function (VACF) on dynamics simulations to get vibrational spectra.

### Harmonic Frequencies

```julia
using YASS

bdys = readSystem("water.xyz")

freqs, vects = getHarmonicFreqs(TIP4Pf(), bdys)
```

`freqs` will be a vector of complex numbers, and `vects` a matrix of floats. Each column in `vects` is the eigenvector associated with a frequency in `freqs`. This means `vects[:,1]` is the eigenvector with frequency `freqs[1]`.

### VACF

The velocity autocorrelation function (VACF) can be used to analyze vibrational modes from molecular dynamics trajectories. `YASS.jl` provides routines to compute the vibrational density of states (VDOS) through the VACF.

To calculate the VACF from a trajectory:

```julia
using YASS

# Read in you water structure
water = readSystem("water.xyz")

# Run MD simulation and get trajectory
traj = run(TIP4Pf(), water, (0.0, 10.0ps), 1fs, NVE())

# Get velocities and masses
vel, mas = getVelMas(traj)

# Create VACF inputs
inp = vacfInps(
    vel,          # Velocity data
    mas,          # Masses
    1e15,         # Sampling frequency (1/fs = 1e15 Hz)
    true,         # Normalize
    Hann,         # Window function
    4,            # Padding factor
    true          # Mirror the data
)

# Calculate VDOS
out = VDOS(inp)
```

The output `out` contains:
- `out.c`: Raw VACF
- `out.C`: Windowed VACF 
- `out.v`: Frequency
- `out.I`: VDOS intensity

You can customize the analysis by:

- Using different window functions (`Hann`, `Welch`, or `HannM`)
- Adjusting the padding factor for FFT
- Enabling/disabling mirroring of the data
- Selecting specific atoms to include using the `atms` keyword argument

For example, to analyze only oxygen atoms:

```julia 
# Get indices of oxygen atoms
O_indices = findall(x -> x == "O", traj.symbols)

# Calculate VDOS for oxygen atoms only
out = VDOS(inp, atms=O_indices)
```

The frequency axis `out.v` is in wavenumbers (cm$^{-1}$) by default.

You can plot the VDOS spectrum using your preferred plotting package:

```julia
using Plots

plot(out.v, out.I, 
     xlabel="Wavenumber (cm-1)", 
     ylabel="Intensity",
     label="VDOS")
```