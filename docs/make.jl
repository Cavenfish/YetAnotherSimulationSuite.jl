using Documenter

makedocs(;
  sitename = "JMD.jl",
  authors = "Brian C. Ferrari",
  pages = [
    "index.md",
    "bodies.md",
    "optimizations.md"
  ]
)

deploydocs(
  repo = "github.com/Cavenfish/JMD.git"
)