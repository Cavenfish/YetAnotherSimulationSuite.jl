# Vibrational Analysis

`JMD.jl` can calculate harmonic vibrational frequencies and eigenvectors. It is also possible to use a velocity autocorrelation function (VACF) on dynamics simulations to get vibrational spectra.

### Harmonic Frequencies

```julia
using JMD

bdys = readSystem("water.xyz")

freqs, vects = getHarmonicFreqs(TIP4Pf(), bdys)
```