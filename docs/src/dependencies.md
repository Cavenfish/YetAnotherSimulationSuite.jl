# Dependencies

YASS relies on several specialized external packages. These packages actively maintained and well trusted within the Julia ecosystem. If this changes, YASS will remove these dependencies and if necessary implement the specialized code in-house. Here the dependencies are listed with links to their repos to give credit to their work but also to provide transperancy.

Packages within the Julia standard library are listed seperately since they are expected to be maintained as well as the Julia langague itself. For each dependency there is a short description of how it is used in YASS, or why it is considered for removal.

**Julia Standard Library Packages**
  - [TOML](https://github.com/JuliaLang/julia/tree/master/stdlib/TOML):
  - [Libdl](https://github.com/JuliaLang/julia/tree/master/stdlib/Libdl):
  - [Statistics](https://github.com/JuliaLang/julia/tree/master/stdlib/Statistics):
  - [Serialization](https://github.com/JuliaLang/julia/tree/master/stdlib/Serialization):
  - [LinearAlgebra](https://github.com/JuliaLang/julia/tree/master/stdlib/LinearAlgebra):

**External Packages**
  - [FFTW](https://github.com/JuliaMath/FFTW.jl):
  - [Optim](https://github.com/JuliaNLSolvers/Optim.jl):
  - [PyCall](https://github.com/JuliaPy/PyCall.jl):
  - [Chemfiles](https://github.com/chemfiles/chemfiles):
  - [Distances](https://github.com/JuliaStats/Distances.jl):
  - [Clustering](https://github.com/JuliaStats/Clustering.jl):
  - [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl):
  - [Distributions](https://github.com/JuliaStats/Distributions.jl):
  - [KernelDensity](https://github.com/JuliaStats/KernelDensity.jl):
  - [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl):

**Considered for Removal**
  - [JLD2](https://github.com/JuliaIO/JLD2.jl): This package is currently only used to load neural network data for potentials. This functionality can be covered by the `Serialization` package, which can reduce the total dependency count. Note, this is a well maintained package and users are encouraged to use it alongside YASS.
  - [DataFrames](https://github.com/JuliaData/DataFrames.jl): This package does not add any functionality to YASS, but rather enchances user exerpience. However, this can be achieved by users using the package alongside YASS, rather than it being a dependency. Note, this is a well maintained package and users are encouraged to use it alongside YASS.