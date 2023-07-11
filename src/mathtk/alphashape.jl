"""
A Simplex in 3D looks like a pyramid. 
Hence it has 4 points associated with it.

My method:
  - Go through all simplexes and form two groups 
    1) perimeter edges
    2) all edges
  - Minimize based on full volume (sum of simplex volume)
"""

struct Edge 
  p1::Vector{Float64}
  p2::Vector{Float64}
end

#Maybe this will be useful?
struct Simplex
  pts::Vector
  r::Float64
  V::Float64
end

function getSimplexVolume(pts)
  N = length(pts)
  n = length(pts[1])

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

  #Get volume
  V = sqrt(det(CM) / c)
  return V
end

function getAlphaVolume(simplexes)

  V = 0.0
  for pts in simplexes
    V += getSimplexVolume(pts)
  end

  return V
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
    b[i] .= dot(pts[i], pts[i])
  end

  #Solve linear equation Ax=b for column vector x
  x = A \ b

  #Get Simplex circumcenter
  c = sum(pts .* x[1:N])

  #Get Simplex circumradius
  r = norm(pts[1] - c)

  return r
end
