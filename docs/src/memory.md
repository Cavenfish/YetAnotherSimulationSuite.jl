# Memory Considerations

Currently, YASS has a roughly 1 GB memory overhead due to dependencies, buffer allocations, and compilation artifacts. This overhead is typical for Julia packages with similar functionality, but may be significant for users with limited memory resources. This overhead is static and does not scale with system size, so larger simulations will see a smaller relative impact. Future optimizations may reduce this overhead.

A few examples of memory usage for different system sizes are provided below, the systems were those used in the [benchmark](https://cavenfish.github.io/YetAnotherSimulationSuite.jl/benchmarks/) section.

## Geometry Optimization

| Number of Atoms | Iterations | Memory Usage (GB) |
|-----------------|------------|-------------------|
| 100             | 1500       | 1.06              |
| 500             | 1500       | 1.06              |
| 1000            | 1500       | 1.06              |


## Harmonic Frequency Calculation

| Number of Atoms | Memory Usage (GB) |
|-----------------|-------------------|
| 100             |  1.07             |
| 500             |  1.16             |
| 1000            |  1.41             |


## Molecular Dynamics Simulation

| Number of Atoms | Time Steps | Memory Usage (GB) |
|-----------------|------------|-------------------|
| 100             | 1000       |  1.12             |
| 500             | 1000       |  1.15             |
| 1000            | 1000       |  1.19             |
| 5000            | 1000       |  1.55             |