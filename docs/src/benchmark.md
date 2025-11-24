# Benchmark

## Disclaimer

Julia often delivers substantial performance gains over Python for numerical and scientific code because it is JIT‑compiled, type‑stable, and generates native LLVM code, so well‑written Julia can approach C/Fortran speeds. However, that speed comes with trade‑offs: just‑in‑time compilation (and package precompilation) introduces startup latency, and Julia’s compilation artifacts and runtime can consume more memory than lightweight Python interpreters. In practice, Julia is most advantageous for long‑running, compute‑intensive workflows; for short scripts or very memory‑constrained environments you should weigh the startup and memory overheads or use precompilation strategies to mitigate them.

## Benchmark Setup

All files used in the benchmark can be found in the [GitHub repo](https://github.com/Cavenfish/YetAnotherSimulationSuite.jl) under the `benchmarks/` directory. All benchmarks were run using an AMD Ryzen 9 3950X 16-Core Processor with 64 GB of RAM running Fedora 43. 

## Benchmark Results

`YetAnotherSimulationSuite.jl` (YASS) was benchmarked against Atomistic Simulation Environment (ASE) for geometry optimization, harmonic frequency calculations, and molecular dynamics simulations. For both packages, a simple Lennard-Jones potential for gold atoms was used.
The results are shown below.

### Short-Running Simulations

In all cases, YASS outperforms ASE once the computational workload is sufficiently large to overcome YASS's initial compilation overhead. The speedups observed can still increase further for longer-running simulations.

#### Geometry Optimization

Geometry optimizations using the LBFGS algorithm were performed on gold clusters containing 100, 500, 1000, and 5000 atoms. The results are shown in the table below.

| Number of Atoms | Iterations | YASS Time (s) | ASE Time (s) | YASS Speedup |
|-----------------|------------|---------------|--------------| -------------|
| 100             | 1500       | 11.8          | **8.4**      | 0.71x        |
| 500             | 1500       | **29.6**      | 64.0         | 2.16x        |
| 1000            | 1500       | **108.7**     | 218.8        | 2.01x        |

#### Harmonic Frequencies

Harmonic frequency calculations were performed on gold clusters containing 100, 500, and 1000 atoms. The results are shown in the table below.

| Number of Atoms | YASS Time (s) | ASE Time (s) | YASS Speedup |
|-----------------|---------------|--------------| -------------|
| 100             | 11.2          | **3.3**      | 0.29x        |
| 500             | **22.7**      | 89.8         | 3.95x        |
| 1000            | **87.1**      | 490.8        | 5.63x        |

#### Molecular Dynamics

Molecular dynamics simulations were performed on gold clusters containing 100, 500, and 1000 atoms for 1000, 2000 and 5000 time steps. The results are shown in the table below.

| Number of Atoms | Time Steps | YASS Time (s) | ASE Time (s) | YASS Speedup |
|-----------------|------------|---------------|--------------| -------------|
| 100             | 1000       | 15.8          | **4.8**      | 0.30x        |
| 100             | 2000       | 16.1          | **9.3**      | 0.58x        |
| 100             | 5000       | **16.6**      | 22.4         | 1.35x        |
| 500             | 1000       | **21.9**      | 32.3         | 1.47x        |
| 500             | 2000       | **27.7**      | 65.9         | 2.38x        |
| 500             | 5000       | **44.2**      | 160.7        | 3.63x        |
| 1000            | 1000       | **35.6**      | 79.4         | 2.23x        |
| 1000            | 2000       | **55.9**      | 165.3        | 2.95x        |
| 1000            | 5000       | **112.5**     | 406.1        | 3.61x        |