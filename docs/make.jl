using Documenter

makedocs(;
  sitename = "YASS.jl",
  authors = "Brian C. Ferrari",
  format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    canonical="https://Cavenfish.github.io/YASS/"
  ),
  pages = [
    "Home" => "index.md",
    "Bodies" => "bodies.md",
    "Optimization" => "optimizations.md",
    "Vibrations" => "vibrations.md",
    "Molecular Dynamics" => ["md/dynamics.md", "md/thermostats.md"],
    "Potentials" => ["potentials/potList.md", "potentials/customPot.md"]
  ]
)

deploydocs(;
  repo = "github.com/Cavenfish/YASS.git",
  devbranch = "dev",
  push_preview = true
)