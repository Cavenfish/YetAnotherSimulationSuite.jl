
struct OrcaOutput
  E::Float64
  μ::Vector{Float64}
  M::Float64
end

function parseOutput(file)
  engStr = "FINAL SINGLE POINT ENERGY"
  μStr   = "Total Dipole Moment"
  MStr   = "Magnitude (Debye)"

  out    = readlines(file)
  
  i = findlast(e -> occursin(engStr, e), out)
  j = findlast(e -> occursin(μStr  , e), out)
  k = findlast(e -> occursin(MStr  , e), out)

  E = split(out[i])[end] |> (x -> parse(Float64, x))
  μ = split(out[j], ':')[end] |> split |> (x -> parse.(Float64, x))
  M = split(out[k])[end] |> (x -> parse(Float64, x))


  OrcaOutput(E, μ, M)
end