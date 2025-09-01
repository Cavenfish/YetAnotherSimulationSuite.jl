"""
Compute the average number density of atoms in a given set.

# Arguments
- `bdys::Vector{MyAtoms}`: Vector of `MyAtoms` objects representing atoms.

# Returns
- `Float64`: Average number density (atoms per unit volume) within the maximum radial extent from the center.
"""
function density(bdys::Vector{MyAtoms})
  o    = center(bdys)
  rmax = [norm(i.r - o) for i in bdys] |> maximum
  V    = 4/3 * pi * rmax^3
  
  length(bdys) / V
end

"""
Compute the average number density for a set of positions relative to a given center.

# Arguments
- `o::Vector{Float64}`: The center point as a vector.
- `pts::Vector{A}`: Vector of position vectors (any subtype of `AbstractArray`).

# Returns
- `Float64`: Average number density (points per unit volume) within the maximum radial extent from the center.
"""
function density(o::Vector{Float64}, pts::Vector{A}) where A<:AbstractArray
  N    = length(pts)
  rmax = [norm(i - o) for i in pts] |> maximum
  V    = 4/3 * pi * rmax^3
  
  N / V
end

"""
Compute the radial distribution function (RDF) from a list of distances and a density.

# Arguments
- `d`: Array of interatomic distances.
- `ρ`: Number density.
- `kwargs...`: Additional keyword arguments passed to `kde_lscv`.

# Returns
- `(x, y)`: Tuple of radius values `x` and normalized RDF values `y`.
"""
function getRdf(d, ρ; kwargs...)
  # get kde and density correction to y
  k  = kde_lscv(d; kwargs...)
  x  = collect(k.x)
  dx = step(k.x)
  C  = @. 4pi * ρ * dx * x ^2
  y  = k.density ./ C 

  x,y
end

"""
Compute the radial distribution function between molecular centers of mass.

# Arguments
- `bdys::Vector{MyAtoms}`: Vector of `MyAtoms` objects.
- `rmol`: Distance threshold for molecular grouping (default 1.2).
- `kwargs...`: Additional keyword arguments for `getRdf`.

# Returns
- `(x, y)`: Tuple of radius values and RDF.
"""
function rdf(bdys::Vector{MyAtoms}; rmol=1.2, kwargs...)
  o    = center(bdys)
  mols = getMols(bdys, rmol)
  pts  = [CoM(bdys[i]) for i in mols]
  ρ    = density(o, pts)
  n    = length(pts)
  d    = [norm(pts[i] - pts[j]) for i = 1:n for j = i+1:n]

  getRdf(d, ρ; kwargs...)
end

"""
Compute the radial distribution function for all pairs of a given atomic species.

# Arguments
- `bdys::Vector{MyAtoms}`: Vector of `MyAtoms` objects.
- `A::String`: Atomic species label.
- `kwargs...`: Additional keyword arguments for `getRdf`.

# Returns
- `(x, y)`: Tuple of radius values and RDF for the given species.
"""
function rdf(bdys::Vector{MyAtoms}, A::String; kwargs...)
  o   = center(bdys)
  i   = findall(e -> e.s == A, bdys)
  pts = [j.r for j in bdys[i]]
  ρ   = density(o, pts)
  n    = length(pts)
  d    = [norm(pts[i] - pts[j]) for i = 1:n for j = i+1:n]

  getRdf(d, ρ; kwargs...)
end

"""
Compute the radial distribution function for pairs of two specified atomic species.

# Arguments
- `bdys::Vector{MyAtoms}`: Vector of `MyAtoms` objects.
- `P::Pair{String, String}`: Pair of atomic species labels (A, B).
- `kwargs...`: Additional keyword arguments for `getRdf`.

# Returns
- `(x, y)`: Tuple of radius values and RDF for the species pair.
"""
function rdf(bdys::Vector{MyAtoms}, P::Pair{String, String}; kwargs...)
  ρ = density(bdys)

  # Get indicies of selected species
  A = findall(e -> e.s == P[1], bdys)
  B = findall(e -> e.s == P[2], bdys)

  # Get distances
  d = [norm(bdys[i].r - bdys[j].r) for i in A for j in B]

  getRdf(d, ρ; kwargs...)
end

"""
Compute the angular distribution function (ADF) between molecular orientations.

# Arguments
- `bdys::Vector{MyAtoms}`: Vector of `MyAtoms` objects.
- `kwargs...`: Additional keyword arguments for `kde_lscv`.

# Returns
- `(x, y)`: Tuple of angle values (in degrees) and ADF.
"""
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
