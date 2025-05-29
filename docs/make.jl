using Documenter

makedocs(;
  sitename = "JMD.jl",
  authors = "Brian C. Ferrari",
  format = Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    canonical="https://Cavenfish.github.io/JMD.jl/"
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

deploydocs(
  repo = "github.com/Cavenfish/JMD.git",
  push_preview = true
)