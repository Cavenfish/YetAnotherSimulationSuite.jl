using Documenter
using YetAnotherSimulationSuite

makedocs(;
  sitename = "YetAnotherSimulationSuite.jl",
  authors = "Brian C. Ferrari",
  format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    canonical="https://Cavenfish.github.io/YetAnotherSimulationSuite.jl/"
  ),
  pages = [
    "Home" => "index.md",
    "Bodies" => "bodies.md",
    "Optimization" => "optimizations.md",
    "Vibrations" => "vibrations.md",
    "Molecular Dynamics" => ["md/dynamics.md", "md/thermostats.md"],
    "Potentials" => ["potentials/potList.md", "potentials/customPot.md"],
    "API" => "api.md",
    "Benchmark" => "benchmark.md",
    "Memory Considerations" => "memory.md",
    "Credits" => "dependencies.md"
  ]
)

deploydocs(;
  repo = "github.com/Cavenfish/YetAnotherSimulationSuite.jl",
  devbranch = "main",
  push_preview = true
)