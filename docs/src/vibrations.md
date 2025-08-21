# Vibrational Analysis

`YASS.jl` can calculate harmonic vibrational frequencies and eigenvectors. It is also possible to use a velocity autocorrelation function (VACF) on dynamics simulations to get vibrational spectra.

### Harmonic Frequencies

```julia
using YASS

bdys = readSystem("water.xyz")

freqs, vects = getHarmonicFreqs(TIP4Pf(), bdys)
```

`freqs` will be a vector of complex numbers, and `vects` a matrix of floats. Each column in `vects` is the eigenvector associated with a frequency in `freqs`. This means `vects[:,1]1` is the eigenvector with frequency `freqs[1]`.

### VACF