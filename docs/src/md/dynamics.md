# Molecular Dynamics Simulations

`YASS.jl` can perform classical molecular dynamics simulations in the NVE and NVT ensemble. 

```julia
using YASS

bdys = readSystem("water.xyz")

ensemble = NVE()

traj = run(TIP4Pf(), bdys, (0.0, 5ps), 0.1fs, ensemble)
```

```julia
using YASS

bdys = readSystem("water.xyz")

ensemble = Berendsen(100.0, 50fs) |> NVT

traj = run(TIP4Pf(), bdys, (0.0, 5ps), 0.1fs, ensemble)
```