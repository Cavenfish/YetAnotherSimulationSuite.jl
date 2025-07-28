
function _pbc!(
  F::AbstractVector, u::AbstractVector, a::Int64, b::Int64,
  func::FT, L::AbstractMatrix, NC::Vector{Int64}, p::PT; cutoff=20.0
) where {FT, PT}
  E = 0.0
  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        r2     = u[b] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        norm(u[a] - r2) > cutoff && continue

        e,f    = func(u[a], r2, p...)
        E     += e
        F[a] .-= f
        F[b] .+= f
      end
    end
  end
  E
end
