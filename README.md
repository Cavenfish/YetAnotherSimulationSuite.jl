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

[![][ddocs-img]][ddocs-url]
[![][ci-img]][ci-url]
[![][codecov-img]][codecov-url]
[![][aqua-img]][aqua-url]


<img src="https://github.com/Cavenfish/YetAnotherSimulationSuite.jl/blob/docs/src/assets/logo.png" alt="Logo" width=350 >

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

```julia
using Pkg
Pkg.add("https://github.com/Cavenfish/YetAnotherSimulationSuite.jl.git")
```

## Contributing

We welcome contributions! Whether it's:

- 🐛 Bug fixes
- ✨ New features
- 📚 Documentation improvements
- 🧪 Additional test cases

If you find YASS useful or just want to show support, please consider starring the repository!

## Development Status

**YASS is currently in pre-release phase.** 

Here's what we're working on:

### Roadmap
- [x] Core functionality cleanup
- [x] Comprehensive documentation
- [ ] Example scripts and tutorials
- [x] Test coverage >70%
- [x] Unitful.jl integration
- [x] Chemfiles.jl integration

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
@misc{yass2023,
  author = {Brian C. Ferrari},
  title = {YetAnotherSimulationSuite.jl},
  year = {2023},
  publisher = {GitHub},
  url = {https://github.com/Cavenfish/YetAnotherSimulationSuite.jl}
}
```