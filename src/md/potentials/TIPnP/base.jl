function getWaterBonds(bdys::Union{MyCell, Vector{MyAtoms}})
  n     = length(bdys)
  bonds = Tuple[]

  for i = 1:3:n
    push!(bonds, (i, i+1))
    push!(bonds, (i, i+2))
  end

  bonds
end

function getMsiteVars!(
  P::T, u::Vector{A}, w1::V, w2::V
) where {T, V <: Vector{Int64}, A <: AbstractVector}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  # Get r vectors
  @. P.r1o = u[h1] - u[o1]
  @. P.r2o = u[h2] - u[o1]
  @. P.r12 = u[h2] - u[h1]
  @. P.r3o = u[h3] - u[o2]
  @. P.r4o = u[h4] - u[o2]
  @. P.r34 = u[h4] - u[h3]

  # Get M1 stuff
  @. P.rbuf = P.r1o + (0.5 * P.r12)
  γ1        = P.drel / norm(P.rbuf)
  P.m1     .= u[o1] .+ P.drel .* (P.rbuf ./ norm(P.rbuf))

  # Get M2 stuff
  @. P.rbuf = P.r3o + (0.5 * P.r34)
  γ2        = P.drel / norm(P.rbuf)
  P.m2     .= u[o2] .+ P.drel .* (P.rbuf ./ norm(P.rbuf))

  γ1, γ2
end

function spreadMforces!(
  F::AbstractVector, Fd::AbstractVector, rid::AbstractVector,
  P::T, w::Vector{Int64}, γ::Float64
) where {T}

  # I'm reusing some 3D vectors to save on
  # memory allocations.

  P.r12 .= dot(rid, Fd) / dot(rid, rid) .* rid
  @. P.r34 = Fd - P.r12

  @. F[w[1]] += Fd - γ * P.r34
  @. F[w[2]] += 0.5 * γ * P.r34
  @. F[w[3]] += 0.5 * γ * P.r34

end

function _getMforces!(
  F::Vector{Af}, u::Vector{Au}, w1::V, w2::V, P::T
) where {Af <: AbstractVector, Au <: AbstractVector, V <: Vector{Int64}, T}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  γ1, γ2 = getMsiteVars!(P, u, w1, w2)

  @. P.r1o = P.m1 - u[o1]
  @. P.r2o = P.m2 - u[o2]

  # H1 -- M2
  E,f     = _Coulomb(u[h1], P.m2, P.Qh, P.Qm)
  f     .*= P.S[1]
  F[h1] .-= f
  spreadMforces!(F, f, P.r2o, P, w2, γ2)

  # H2 -- M2
  e,f     = _Coulomb(u[h2], P.m2, P.Qh, P.Qm)
  f     .*= P.S[1]
  E      += e
  F[h2] .-= f
  spreadMforces!(F, f, P.r2o, P, w2, γ2)

  # H3 -- M1
  e,f     = _Coulomb(u[h3], P.m1, P.Qh, P.Qm)
  f     .*= P.S[1]
  E      += e
  F[h3] .-= f
  spreadMforces!(F, f, P.r1o, P, w1, γ1)

  # H4 -- M1
  e,f     = _Coulomb(u[h4], P.m1, P.Qh, P.Qm)
  f     .*= P.S[1]
  E      += e
  F[h4] .-= f
  spreadMforces!(F, f, P.r1o, P, w1, γ1)

  # M1 -- M2
  e,f       = _Coulomb(P.m1, P.m2, P.Qm, P.Qm)
  @. P.Fbuf = f * P.S[1]
  E        += e
  spreadMforces!(F, P.Fbuf, P.r2o, P, w2, γ2)
  P.Fbuf .*= -1
  spreadMforces!(F, P.Fbuf, P.r1o, P, w1, γ1)

  E
end

function tip4pf_pbc!(
  F::Vector{Af}, u::Vector{Au}, w1::Vi, w2::Vi,
  NC::Vector{Int64}, L::AbstractMatrix, P::T
) where {Af, Au, Vi <: Vector{Int64}, T}

  E = 0.0
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        @. P.rbuf2 = (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        @. P.o2t  = u[o2] + P.rbuf2
        doo       = norm(P.o2t - u[o1])
        @. P.rbuf = (P.o2t - u[o1]) / doo
        switchLR!(P.S, doo, P.rs, P.rc)

        if iszero(P.S)
          continue
        end

        @. P.h3t = u[h3] + P.rbuf2
        @. P.h4t = u[h4] + P.rbuf2

        # O1-O2 VDW
        e,f      = _vdw(u[o1], P.o2t, P.ϵoo, P.σoo)
        E        += e * P.S[1]
        P.Fbuf   .= -P.S[2] * e .* P.rbuf
        @. F[o1] -= f * P.S[1]
        @. F[o2] += f * P.S[1]

        # H1-H3 Coulomb
        e,f       = _Coulomb(u[h1], P.h3t, P.Qh, P.Qh)
        E        += e * P.S[1]
        P.Fbuf  .+= -P.S[2] * e .* P.rbuf
        @. F[h1] -= f * P.S[1]
        @. F[h3] += f * P.S[1]

        # H2-H4 Coulomb
        e,f       = _Coulomb(u[h1], P.h4t, P.Qh, P.Qh)
        E        += e * P.S[1]
        P.Fbuf  .+= -P.S[2] * e .* P.rbuf
        @. F[h1] -= f * P.S[1]
        @. F[h4] += f * P.S[1]

        # H2-H3 Coulomb
        e,f       = _Coulomb(u[h2], P.h3t, P.Qh, P.Qh)
        E        += e * P.S[1]
        P.Fbuf  .+= -P.S[2] * e .* P.rbuf
        @. F[h2] -= f * P.S[1]
        @. F[h3] += f * P.S[1]

        # H2-H4 Coulomb
        e,f       = _Coulomb(u[h2], P.h4t, P.Qh, P.Qh)
        E        += e * P.S[1]
        P.Fbuf  .+= -P.S[2] * e .* P.rbuf
        @. F[h2] -= f * P.S[1]
        @. F[h4] += f * P.S[1]

        # Spread dS force
        F[o1] .-= P.Fbuf
        F[o2] .+= P.Fbuf
      end
    end
  end      

  E
end

function pbc_Mforces!(
  F::AbstractVector, u::AbstractVector, w1::V, w2::V, 
  NC::V, L::AbstractMatrix, P::T
) where {V <: Vector{Int64}, T}
  
  E = 0.0
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  γ1, γ2 = getMsiteVars!(P, u, w1, w2)

  @. P.r1o = P.m1 - u[o1]
  @. P.r2o = P.m2 - u[o2]

  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        @. P.rbuf2 = (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        @. P.o2t  = u[o2] + P.rbuf2
        doo       = norm(P.o2t - u[o1])
        @. P.rbuf = (P.o2t - u[o1]) / doo
        switchLR!(P.S, doo, P.rs, P.rc)

        if iszero(P.S)
          continue
        end

        @. P.m2t = P.m2  + P.rbuf2
        @. P.h3t = u[h3] + P.rbuf2
        @. P.h4t = u[h4] + P.rbuf2

        # H1 -- M2
        e,f     = _Coulomb(u[h1], P.m2t, P.Qh, P.Qm)
        E      += e * P.S[1]
        P.Fbuf .= -P.S[2] * e .* P.rbuf
        f     .*= P.S[1]
        F[h1] .-= f
        spreadMforces!(F, f, P.r2o, P, w2, γ2)

        # H2 -- M2
        e,f      = _Coulomb(u[h2], P.m2t, P.Qh, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        F[h2]  .-= f
        spreadMforces!(F, f, P.r2o, P, w2, γ2)

        # H3 -- M1
        e,f      = _Coulomb(P.h3t, P.m1, P.Qh, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        F[h3]  .-= f
        spreadMforces!(F, f, P.r1o, P, w1, γ1)

        # H4 -- M1
        e,f      = _Coulomb(P.h4t, P.m1, P.Qh, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        F[h4]  .-= f
        spreadMforces!(F, f, P.r1o, P, w1, γ1)

        # M1 -- M2
        e,f      = _Coulomb(P.m1, P.m2t, P.Qm, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        spreadMforces!(F, f, P.r2o, P, w2, γ2)
        f .*= -1
        spreadMforces!(F, f, P.r1o, P, w1, γ1)

        # Spread dS force
        F[o1] .-= P.Fbuf
        F[o2] .+= P.Fbuf
      end
    end
  end

  E
end