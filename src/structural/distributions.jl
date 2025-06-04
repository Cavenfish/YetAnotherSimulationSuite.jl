
function density(bdys::Vector{MyAtoms})
  o    = center(bdys)
  rmax = [norm(i.r - o) for i in bdys] |> maximum
  V    = 4/3 * pi * rmax^3
  
  length(bdys) / V
end

function density(o::Vector{Float64}, pts::Vector{A}) where A<:AbstractArray
  N    = length(pts)
  rmax = [norm(i - o) for i in pts] |> maximum
  V    = 4/3 * pi * rmax^3
  
  N / V
end

function getRdf(d, ρ; kwargs...)
  # get kde and density correction to y
  k  = kde_lscv(d; kwargs...)
  x  = collect(k.x)
  dx = step(k.x)
  C  = @. 4pi * ρ * dx * x ^2
  y  = k.density ./ C 

  x,y
end

function rdf(bdys::Vector{MyAtoms}; rmol=1.2, kwargs...)
  o    = center(bdys)
  mols = getMols(bdys, rmol)
  pts  = [CoM(bdys[i]) for i in mols]
  ρ    = density(o, pts)
  n    = length(pts)
  d    = [norm(pts[i] - pts[j]) for i = 1:n for j = i+1:n]

  getRdf(d, ρ; kwargs...)
end

function rdf(bdys::Vector{MyAtoms}, A::String; kwargs...)
  o   = center(bdys)
  i   = findall(e -> e.s == A, bdys)
  pts = [j.r for j in bdys[i]]
  ρ   = density(o, pts)
  n    = length(pts)
  d    = [norm(pts[i] - pts[j]) for i = 1:n for j = i+1:n]

  getRdf(d, ρ; kwargs...)
end

function rdf(bdys::Vector{MyAtoms}, P::Pair{String, String}; kwargs...)
  ρ = density(bdys)

  # Get indicies of selected species
  A = findall(e -> e.s == P[1], bdys)
  B = findall(e -> e.s == P[2], bdys)

  # Get distances
  d = [norm(bdys[i].r - bdys[j].r) for i in A for j in B]

  getRdf(d, ρ; kwargs...)
end

function adf(bdys::Vector{MyAtoms}; kwargs...)
  ρ = density(bdys)

  m = getMols(bdys)
  n = length(m)
  θ = [getAngleCO(bdys[m[i]], bdys[m[j]]) for i = 1:n for j = i+1:n] .* (360/2pi)
  k = kde_lscv(θ; kwargs...)
  x = collect(k.x)
  y = k.density
  return x,y
end
