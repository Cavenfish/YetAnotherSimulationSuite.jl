"""
TIP4P/2005f 

Based on:
González, M. A., & Abascal, J. L. (2011). A flexible model for water based on
TIP4P/2005. The Journal of chemical physics, 135(22).

Link:
https://pubs.aip.org/aip/jcp/article/135/22/224516/190786
"""

TIP4Pf(; constraints=nothing) = Calculator(TIP4Pf; EF=TIP4Pf!, constraints=constraints)

struct _TIP4P_PotVars{
  F<:Float64, AV3D<:AbstractVector, AV2D<:AbstractVector
} <: PotVars
  D::F
  a::F
  req::F
  K::F
  θeq::F
  ϵoo::F
  σoo::F
  drel::F
  Qh::F
  Qm::F
  rs::F
  rc::F
  Fbuf::AV3D
  rbuf::AV3D
  rbuf2::AV3D
  m1::AV3D
  m2::AV3D
  S::AV2D
  o2t::AV3D
  h3t::AV3D
  h4t::AV3D
  m2t::AV3D
  r1o::AV3D
  r2o::AV3D
  r12::AV3D
  r3o::AV3D
  r4o::AV3D
  r34::AV3D
end

#PotVar building function 
TIP4Pf(bdys::Union{MyCell, Vector{MyAtoms}}) = _TIP4P_PotVars(
  4.48339,    # eV
  2.287,      # \AA
  0.9419,     # \AA
  3.81209321, # eV rad^-2  1.16123e-3 * (2pi/360)^2
  1.87448,    # rad        107.4  * (2pi/360)
  8.03e-3,    # eV
  3.1644,     # \AA
  0.13194,    # \AA
  2.1113635,  # 
  -2 * 2.1113635, #
  11.0,       # \AA
  12.0,       # \AA
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{2}(zeros(2)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
)

function TIP4Pf!(F, u, p)
  E = 0.0
  P = p.potVars

  for mol in p.mols
    o, h1, h2 = mol

    E += _Morse!(F, u, o, h1, P.D, P.a, P.req)
    E += _Morse!(F, u, o, h2, P.D, P.a, P.req)
    E += _harmonicBondAngle!(F, u, h1, o, h2, P.K, P.θeq)
  end

  for i = 1:length(p.mols)
    o1, h1, h2 = p.mols[i]

    for j = i+1:length(p.mols)
      o2, h3, h4 = p.mols[j]

      @. P.rbuf2 = u[o2] - u[o1]
      doo = norm(P.rbuf2)
      switchLR!(P.S, doo, P.rs, P.rc)

      if iszero(P.S)
        continue
      end

      e = 0.0

      # O-O vdw Interaction
      e += _vdw!(F, u, P.Fbuf, P.rbuf, o1, o2, P.ϵoo, P.σoo; S=P.S[1])

      # H-H Coulomb Interactions
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h1, h3, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h1, h4, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h2, h3, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h2, h4, P.Qh, P.Qh; S=P.S[1])

      e += _getMforces!(F, u, p.mols[i], p.mols[j], P)

      E += P.S[1] * e

      if P.S[2] != 0.0
        @. P.rbuf = P.rbuf2 / doo
        P.Fbuf .= -P.S[2] * e .* P.rbuf
        F[o1] .-= P.Fbuf
        F[o2] .+= P.Fbuf
      end
    end
  end

  if any(p.PBC)
    NC    = p.NC .* p.PBC
    lat   = isa(p, optVars) ? p.lattice : p.ensemble.lattice

    for i = 1:length(p.mols)
      a = p.mols[i]

      for j = i+1:length(p.mols)
        b = p.mols[j]

        # vdw and H-H Coulomb
        E += tip4pf_pbc!(F, u, a, b, NC, lat, P)

        # M-site Coulomb
        E += pbc_Mforces!(F, u, a, b, NC, lat, P)
      end
    end
  end

  E
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

  # Get a
  a1 = 1 / (1 + (norm(P.r2o) / norm(P.r1o)))
  a2 = 1 / (1 + (norm(P.r4o) / norm(P.r3o)))

  # Get M1 stuff
  @. P.rbuf = P.r1o + (a1 * P.r12)
  γ1        = P.drel / norm(P.rbuf)
  P.m1     .= u[o1] .+ P.drel .* (P.rbuf ./ norm(P.rbuf))

  # Get M2 stuff
  @. P.rbuf = P.r3o + (a2 * P.r34)
  γ2        = P.drel / norm(P.rbuf)
  P.m2     .= u[o2] .+ P.drel .* (P.rbuf ./ norm(P.rbuf))

  a1, a2, γ1, γ2
end

function spreadMforces!(
  F::AbstractVector, Fd::AbstractVector, rid::AbstractVector,
  P::T, w::Vector{Int64}, γ::Float64, a::Float64
) where {T}

  # I'm reusing some 3D vectors to save on
  # memory allocations.

  P.r12 .= dot(rid, Fd) / dot(rid, rid) .* rid
  @. P.r34 = Fd - P.r12

  @. F[w[1]] += Fd - γ * P.r34
  @. F[w[2]] += (1 - a) * γ * P.r34
  @. F[w[3]] += a * γ * P.r34

end

function _getMforces!(
  F::Vector{Af}, u::Vector{Au}, w1::V, w2::V, P::T
) where {Af <: AbstractVector, Au <: AbstractVector, V <: Vector{Int64}, T}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  a1, a2, γ1, γ2 = getMsiteVars!(P, u, w1, w2)

  @. P.r1o = P.m1 - u[o1]
  @. P.r2o = P.m2 - u[o2]

  # H1 -- M2
  E,f     = _Coulomb(u[h1], P.m2, P.Qh, P.Qm)
  f     .*= P.S[1]
  F[h1] .-= f
  spreadMforces!(F, f, P.r2o, P, w2, γ2, a2)

  # H2 -- M2
  e,f     = _Coulomb(u[h2], P.m2, P.Qh, P.Qm)
  f     .*= P.S[1]
  E      += e
  F[h2] .-= f
  spreadMforces!(F, f, P.r2o, P, w2, γ2, a2)

  # H3 -- M1
  e,f     = _Coulomb(u[h3], P.m1, P.Qh, P.Qm)
  f     .*= P.S[1]
  E      += e
  F[h3] .-= f
  spreadMforces!(F, f, P.r1o, P, w1, γ1, a1)

  # H4 -- M1
  e,f     = _Coulomb(u[h4], P.m1, P.Qh, P.Qm)
  f     .*= P.S[1]
  E      += e
  F[h4] .-= f
  spreadMforces!(F, f, P.r1o, P, w1, γ1, a1)

  # M1 -- M2
  e,f       = _Coulomb(P.m1, P.m2, P.Qm, P.Qm)
  @. P.Fbuf = f * P.S[1]
  E        += e
  spreadMforces!(F, P.Fbuf, P.r2o, P, w2, γ2, a2)
  P.Fbuf .*= -1
  spreadMforces!(F, P.Fbuf, P.r1o, P, w1, γ1, a1)

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

  a1, a2, γ1, γ2 = getMsiteVars!(P, u, w1, w2)

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
        spreadMforces!(F, f, P.r2o, P, w2, γ2, a2)

        # H2 -- M2
        e,f      = _Coulomb(u[h2], P.m2t, P.Qh, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        F[h2]  .-= f
        spreadMforces!(F, f, P.r2o, P, w2, γ2, a2)

        # H3 -- M1
        e,f      = _Coulomb(P.h3t, P.m1, P.Qh, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        F[h3]  .-= f
        spreadMforces!(F, f, P.r1o, P, w1, γ1, a1)

        # H4 -- M1
        e,f      = _Coulomb(P.h4t, P.m1, P.Qh, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        F[h4]  .-= f
        spreadMforces!(F, f, P.r1o, P, w1, γ1, a1)

        # M1 -- M2
        e,f      = _Coulomb(P.m1, P.m2t, P.Qm, P.Qm)
        E       += e * P.S[1]
        P.Fbuf .+= -P.S[2] * e .* P.rbuf
        f      .*= P.S[1]
        spreadMforces!(F, f, P.r2o, P, w2, γ2, a2)
        f .*= -1
        spreadMforces!(F, f, P.r1o, P, w1, γ1, a1)

        # Spread dS force
        F[o1] .-= P.Fbuf
        F[o2] .+= P.Fbuf
      end
    end
  end

  E
end