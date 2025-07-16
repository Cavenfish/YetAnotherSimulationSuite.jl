
function _pbc!(F, u, a, b, func, L, NC, p; cutoff=20.0)
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
      end
    end
  end
  E
end

function pbc_vdw!(
  F::Vector{Vf}, u::Vector{Vu}, i::Int64, js::Vector{Int64}, 
  ϵij::Float64, σij::Float64, NC::Vector{Int64}, L::AbstractMatrix;
  cutoff=20.0
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  E = 0.0

  for j in js
    E += _pbc!(F, u, i, j, _vdw, L, NC, (ϵij, σij); cutoff=cutoff)
  end

  E
end

function pbc_Coulomb!(
  F::Vector{Vf}, u::Vector{Vu}, i::Int64, js::Vector{Int64}, 
  Qi::Float64, Qj::Float64, NC::Vector{Int64}, L::AbstractMatrix;
  cutoff=20.0
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  E = 0.0

  for j in js
    E += _pbc!(F, u, i, j, _Coulomb, L, NC, (Qi, Qj); cutoff=cutoff)
  end

  E
end
