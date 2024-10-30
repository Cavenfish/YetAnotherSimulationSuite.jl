"""
A Simplex in 3D looks like a pyramid. 
Hence it has 4 points associated with it.

The edge of a 3D Simplex is a 2D triangle.

The volume of an edge is really an area

My method:
  - Go through all simplexes and form two groups 
    1) perimeter edges
    2) all edges
  - Minimize based on full volume (sum of simplex volume)
"""

struct Edge 
  pts::Vector{Vector{Float64}}
  indices::Vector{Int64}
  V::Float64
end

struct Simplex
  pts::Vector{Vector{Float64}}
  indices::Vector{Int64}
  r::Float64
  V::Float64
end

struct AlphaShape
  pts::Vector
  perimeter::Vector
  simplexes::Vector
  area::Float64
end

function getEdges(simp)
  n = length(simp.indices)
  e = Edge[]

  for i in 1:n
    for j in i+1:n
      for k in j+1:n
        pts = [simp.pts[q] for q in [i,j,k]]
        ind = [simp.indices[q] for q in [i,j,k]]
        V   = getSimplexVolume(pts)
        push!(e, Edge(pts, ind, V))
      end
    end
  end

  e
end

function getSimplexes(pts)
  f  = "qhull d Qt Qbb Qc Qz Q12"
  DT = delaunay(pts, f)
  n  = size(DT)[2]

  simplexes = Simplex[]
  for i in 1:n
    p = pts[DT[:,i]]
    j = DT[:,i]
    r = getSimplexRadius(p)
    V = getSimplexVolume(p)
    push!(simplexes, Simplex(p,j,r,V))
  end

  simplexes
end

function getSimplexVolume(pts)
  N = length(pts)
  n = N - 1

  #First we build the CayleyMenger Matrix
  CM      = ones(N+1, N+1)
  CM[1,1] = 0.0
  for i in 2:N+1
    for j in 2:N+1
      CM[i,j] = norm(pts[i-1] - pts[j-1])^2
    end
  end

  #Get constant
  c = (-1)^(n+1) * factorial(n)^2 * 2^n

  sqrt(det(CM) / c)
end

function getSimplexRadius(pts)
  N = length(pts)
  n = length(pts[1])

  V = zeros(n,N)
  for i in 1:N
    V[:,i] .= pts[i]
  end

  #Build A matrix
  A            = ones(N+1, N+1)
  A[N,N]       = 0.0 
  A[1:N, 1:N] .= 2 * transpose(V) * V

  #Build b column vector
  b = ones(5)
  for i in 1:N
    b[i] = dot(pts[i], pts[i])
  end

  #Solve linear equation Ax=b for column vector x
  x = A \ b

  #Get Simplex circumcenter
  c = sum(pts .* x[1:N])

  norm(pts[1] - c)
end

function getAlpha(pts)
  function f(x)
    A = alphashape(pts; α=x[1])
    V = sum([i.V for i in A.simplexes])
    p = [j for i in A.simplexes for j in i.pts]
    
    l = setdiff(pts, p) |> length
    n = length(A.perimeter)
    m = vcat(A.perimeter...) |> unique |> length
    r = (A.area / V)*(n/m) + (l / m) * 1e65
    
    r
  end

  res = optimize(f, 0.0, 1.0)
  
  res.minimizer[1]
end

function alphashape(pts; α=nothing)

  if α == nothing
    α = getAlpha(pts)
  end
  
  area      = 0.0
  allEdges  = []
  perimeter = []
  simplexes = getSimplexes(pts)

  for simp in simplexes

    #Filter test
    if simp.r < 1/α

      #pull all edges from simplex
      edges = getEdges(simp)

      for edge in edges

        ind = sort(edge.indices)

        # Perimeter edges can only exist once,
        # so here we check if the edge has been
        # seen before. If yes, remove it from 
        # perimeter list (only if it is there)
        if ind in allEdges
          i  = findfirst(e -> e == ind, perimeter)
          
          if i != nothing
            popat!(perimeter, i)
            area -= edge.V
          end

        else
          push!(perimeter, ind)
          area += edge.V
        end

        push!(allEdges, ind)
      end
    else
      i = findfirst(e -> e == simp, simplexes)
      popat!(simplexes, i)
    end
  end

  AlphaShape(pts, perimeter, simplexes, area)
end