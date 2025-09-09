"""
SPC-F

Based on:
Toukan, Kahled, and Aneesur Rahman. "Molecular-dynamics study of atomic
motions in water." Physical Review B 31.5 (1985): 2643.

Link:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.31.2643
"""
SPCF(; constraints=nothing) = Calculator(SPCF; EF=SPCF!, constraints=constraints)

struct _SPCF_PotVars{
  F<:Float64, SV3D<:AbstractVector, AV3D<:AbstractVector, AV2D<:AbstractVector
} <: PotVars
  Kb::F
  req::F
  Kθ::F
  θeq::F
  σ::F
  ϵ::F
  Qo::F
  Qh::F
  rs::F
  rc::F
  Fbuf::AV3D
  rbuf::AV3D
  rbuf2::AV3D
  S::AV2D
  o2t::AV3D
  h3t::AV3D
  h4t::AV3D
end

SPCF(bdys::Union{Vector{MyAtoms}, MyCell}) = _SPCF_PotVars(
  48.05913,
  1.0,
  3.97,
  1.910611,
  3.145, 
  0.007,
  -2.959855,
  0.5 * 2.959855,
  9.0,
  10.0,
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{2}(zeros(2)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
)

function SPCF!(F, u, p)
  E   = 0.0
  P   = p.potVars
  NC  = p.NC .* p.PBC
  lat = isa(p, optVars) ? p.lattice : p.ensemble.lattice

  for mol in p.mols
    o, h1, h2 = mol

    E += _harmonicBond!(F, u, P.Fbuf, P.rbuf, o, h1, P.Kb, P.req)
    E += _harmonicBond!(F, u, P.Fbuf, P.rbuf, o, h2, P.Kb, P.req)
    E += _harmonicBondAngle!(F, u, h1, o, h2, P.Kθ, P.θeq)
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

      e  = 0.0

      # O-O Interactions
      e += _vdw!(F, u, P.Fbuf, P.rbuf, o1, o2, P.ϵ, P.σ; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, o1, o2, P.Qo, P.Qo; S=P.S[1])

      # H-H Coulomb
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h1, h3, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h1, h4, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h2, h3, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, h2, h4, P.Qh, P.Qh; S=P.S[1])

      # H-O Coulomb
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, o1, h3, P.Qo, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, o1, h4, P.Qo, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, o2, h1, P.Qo, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, o2, h2, P.Qo, P.Qh; S=P.S[1])

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
    for i = 1:length(p.mols)
      for j = i+1:length(p.mols)
        E += spcf_pbc!(F, u, p.mols[i], p.mols[j], NC, lat, P)
      end
    end
  end

  E
end

function spcf_pbc!(
  F::Vector{Af}, u::Vector{Au}, w1::Vi, w2::Vi,
  NC::Vector{Int64}, L::AbstractMatrix, p::P
) where {Af, Au, Vi <: Vector{Int64}, P}

  E = 0.0
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        @. p.rbuf2 = (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        @. p.o2t = u[o2] + p.rbuf2
        roo = norm(u[o1] - p.o2t)
        p.rbuf .= (p.o2t - u[o1]) / roo
        switchLR!(p.S, roo, p.rs, p.rc)

        if iszero(p.S)
          continue
        end

        @. p.h3t = u[h3] + p.rbuf2
        @. p.h4t = u[h4] + p.rbuf2

        # O-O VDW
        e,f    = _vdw(u[o1], p.o2t, p.ϵ, p.σ)
        E     += e * p.S[1]
        p.Fbuf   .= -p.S[2] * e * p.rbuf
        @. F[o1] -= f * p.S[1]
        @. F[o2] += f * p.S[1]

        # O-O Coulomb
        e,f    = _Coulomb(u[o1], p.o2t, p.Qo, p.Qo)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[o1] -= f * p.S[1]
        @. F[o2] += f * p.S[1]

        # H1-H3 Coulomb
        e,f    = _Coulomb(u[h1], p.h3t, p.Qh, p.Qh)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[h1] -= f * p.S[1]
        @. F[h3] += f * p.S[1]

        # H1-H4 Coulomb
        e,f    = _Coulomb(u[h1], p.h4t, p.Qh, p.Qh)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[h1] -= f * p.S[1]
        @. F[h4] += f * p.S[1]

        # H2-H3 Coulomb
        e,f    = _Coulomb(u[h2], p.h3t, p.Qh, p.Qh)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[h2] -= f * p.S[1]
        @. F[h3] += f * p.S[1]

        # H2-H4 Coulomb
        e,f    = _Coulomb(u[h2], p.h4t, p.Qh, p.Qh)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[h2] -= f * p.S[1]
        @. F[h4] += f * p.S[1]

        # H1-O2 Coulomb
        e,f    = _Coulomb(u[h1], p.o2t, p.Qh, p.Qo)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[h1] -= f * p.S[1]
        @. F[o2] += f * p.S[1]

        # H2-O2 Coulomb
        e,f    = _Coulomb(u[h2], p.o2t, p.Qh, p.Qo)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[h2] -= f * p.S[1]
        @. F[o2] += f * p.S[1]

        # H3-O1 Coulomb
        e,f    = _Coulomb(u[o1], p.h3t, p.Qo, p.Qh)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[o1] -= f * p.S[1]
        @. F[h3] += f * p.S[1]

        # H4-O1 Coulomb
        e,f    = _Coulomb(u[o1], p.h4t, p.Qo, p.Qh)
        E      += e * p.S[1]
        p.Fbuf    .+= -p.S[2] * e * p.rbuf
        @. F[o1] -= f * p.S[1]
        @. F[h4] += f * p.S[1]

        # Spread dS force
        F[o1] .-= p.Fbuf
        F[o2] .+= p.Fbuf
      end
    end
  end

  E
end