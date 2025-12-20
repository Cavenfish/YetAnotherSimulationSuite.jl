[ci-img]: https://github.com/Cavenfish/YetAnotherSimulationSuite.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/Cavenfish/YetAnotherSimulationSuite.jl/actions/workflows/CI.yml

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[codecov-img]: https://codecov.io/github/Cavenfish/YetAnotherSimulationSuite.jl/branch/main/graph/badge.svg
[codecov-url]: https://app.codecov.io/github/Cavenfish/YetAnotherSimulationSuite.jl

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://cavenfish.github.io/YetAnotherSimulationSuite.jl/stable/

[ddocs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[ddocs-url]: https://cavenfish.github.io/YetAnotherSimulationSuite.jl/dev/

[joss-img]: https://joss.theoj.org/papers/10.21105/joss.09480/status.svg
[joss-url]: https://doi.org/10.21105/joss.09480

[![][docs-img]][docs-url]
[![][ddocs-img]][ddocs-url]
[![][ci-img]][ci-url]
[![][codecov-img]][codecov-url]
[![][aqua-img]][aqua-url]
[![][joss-img]][joss-url]

<img src="https://github.com/Cavenfish/YetAnotherSimulationSuite.jl/blob/main/docs/src/assets/logo.png" alt="Logo" width=350 >

# YetAnotherSimulationSuite.jl (YASS)

YASS is a modern, flexible atomic simulation suite written in Julia. It aims to provide:

- 🎯 Simple, and intuitive API
- ⚡ High performance native Julia implementation
- 🔧 Easy extensibility for custom methods
- 📦 Built-in potentials and analysis tools

## Quick Start

```julia
using YetAnotherSimulationSuite

# Read molecule
water = readSystem("water.xyz")

# Run 10ps NVE simulation
traj = run(TIP4Pf(), water, (0.0, 10.0ps), 1.0fs, NVE())

# Analyze results
energies = [img.energy for img in traj.images]
temperatures = [img.temp for img in traj.images]
```

## Features

- 🧪 Multiple molecular dynamics ensembles (NVE, NVT)
- 🔬 Built-in analysis tools (RDF, VACF)
- ⚛️ Geometry and cell optimizations
- 📊 Common water models (TIP4P/2005f, SPC-F) 
- 💻 Easy-to-extend architecture
- 🚄 High performance through Julia's native speed
- 📝 Comprehensive documentation

## Installation

```
pkg> add YetAnotherSimulationSuite
```

## Performance

Julia often delivers substantial performance gains over Python for numerical and scientific code because it is JIT‑compiled, type‑stable, and generates native LLVM code, so well‑written Julia can approach C/Fortran speeds. However, that speed comes with trade‑offs: just‑in‑time compilation (and package precompilation) introduces startup latency, and Julia’s compilation artifacts and runtime can consume more memory than lightweight Python interpreters. In practice, Julia is most advantageous for long‑running, compute‑intensive workflows; for short scripts or very memory‑constrained environments you should weigh the startup and memory overheads or use precompilation strategies to mitigate them.

A benchmark comparing YASS to other similar packages can be found in the [documentation](https://cavenfish.github.io/YetAnotherSimulationSuite.jl/stable/benchmark/).

## Memory Considerations

Currently, YASS has a roughly 1 GB memory overhead due to dependencies, buffer allocations, and compilation artifacts. This overhead is typical for Julia packages with similar functionality, but may be significant for users with limited memory resources. This overhead is static and does not scale with system size, so larger simulations will see a smaller relative impact. Future optimizations may reduce this overhead.

A few examples of memory usage for different system sizes can be found in the [documentation](https://cavenfish.github.io/YetAnotherSimulationSuite.jl/stable/memory/).

## Contributing

We welcome contributions! Whether it's:

- 🐛 Bug fixes
- ✨ New features
- 📚 Documentation improvements
- 🧪 Additional test cases

If you find YASS useful or just want to show support, please consider starring the repository!

## Development Status

Here's what we're working on:

### Upcoming Features
- [ ] NPT ensemble simulations
- [ ] Anharmonic vibrational analysis
- [ ] Path integral molecular dynamics (PIMD)
- [ ] Additional analysis tools

## License

YASS is MIT licensed. See [LICENSE](LICENSE) for details.

## Citation

If you use YASS in your research, please cite:

```bibtex
@article{Ferrari2025,
  doi = {10.21105/joss.09480},
  url = {https://doi.org/10.21105/joss.09480},
  year = {2025},
  publisher = {The Open Journal},
  volume = {10},
  number = {116},
  pages = {9480},
  author = {Ferrari, Brian C.},
  title = {YetAnotherSimulationSuite.jl: An Atomic Simulation Suite in Julia},
  journal = {Journal of Open Source Software}
}
```