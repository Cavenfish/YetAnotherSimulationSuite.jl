[ci-img]: https://github.com/Cavenfish/YASS.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/Cavenfish/YASS.jl/actions/workflows/CI.yml

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[codecov-img]: https://codecov.io/github/Cavenfish/YASS.jl/branch/main/graph/badge.svg
[codecov-url]: https://app.codecov.io/github/Cavenfish/YASS.jl

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://cavenfish.github.io/YASS.jl/stable/

[ddocs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[ddocs-url]: https://cavenfish.github.io/YASS.jl/dev/

# Yet Another Simulation Suite (YASS.jl)

[![][ddocs-img]][ddocs-url]
[![][ci-img]][ci-url]
[![][codecov-img]][codecov-url]
[![][aqua-img]][aqua-url]


**`YASS.jl` is still in the pre-release phase.** The roadmap below tracks the development process of `YASS.jl`.

`YASS.jl` aims to offer users a simple, intuitive and easy-to-use molecular dynamics enviornment. It draws inspiration from Python's [ASE](https://wiki.fysik.dtu.dk/ase/index.html), but is intended to be faster and offer users more flexibility. The flexibility comes from the relative ease with which users can add their own methods to dynamics or other components of `YASS.jl`.

# Roadmap to Becoming a True Package

  - [x] Clear out niche code
  - [ ] Make documentation
  - [ ] Add example scripts
  - [x] Get codecov above 70%
  - [ ] Rewrite code to use Unitful.jl?
  - [x] Use Chemfiles for file IO

# Features to Add

  - [ ] Add NPT simulations
  - [ ] Add anharmonic vibrational analysis
  - [ ] Add PIMD methods