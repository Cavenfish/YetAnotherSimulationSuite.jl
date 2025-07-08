
# Maybe????
macro pbc(func)
  
end

function _pbc!(F, u, a, b, func, L, NC, p)
  E = 0.0
  for i = 1:3
    for j = -NC[i]:NC[i]
      j == 0 && continue
      
      r2    = u[b] + (L[i, :] * j)
      e,f   = func(u[a], r2, p...)
      E    += e
      F[a] -= f
    end
  end
  E
end

function pbc_vdw!(
  F::Vector{Vf}, u::Vector{Vu}, i::Int64, js::Vector{Int64}, 
  eij::Float64, oij::Float64, NC::Vector{Int64}, L::AbstractMatrix
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  E = 0.0

  for j in js
    E += _pbc!(F, u, i, j, _vdw, L, NC, (eij, oij))
  end

  E
end

function pbc_Coulomb!(
  F::Vector{Vf}, u::Vector{Vu}, i::Int64, js::Vector{Int64}, 
  Qi::Float64, Qj::Float64, NC::Vector{Int64}, L::AbstractMatrix
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  E = 0.0

  for j in js
    E += _pbc!(F, u, i, j, _Coulomb, L, NC, (Qi, Qj))
  end

  E
end
