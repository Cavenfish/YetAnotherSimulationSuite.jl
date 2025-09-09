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
  F<:Float64, AV3D<:AbstractVector, AV2D<:AbstractVector
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
  lat = isa(p, optVars) ? p.lattice : p.ensemble.lattice

  for mol in p.mols
    o, h1, h2 = mol

    E += _harmonicBond!(F, u, P.Fbuf, P.rbuf, lat, o, h1, P.Kb, P.req)
    E += _harmonicBond!(F, u, P.Fbuf, P.rbuf, lat, o, h2, P.Kb, P.req)
    E += _harmonicBondAngle!(F, u, P.rbuf, P.rbuf2, lat, h1, o, h2, P.Kθ, P.θeq)
  end

  for i = 1:length(p.mols)
    o1, h1, h2 = p.mols[i]

    for j = i+1:length(p.mols)
      o2, h3, h4 = p.mols[j]

      doo = pbcVec!(P.rbuf2, u[o1], u[o2], lat)
      switchLR!(P.S, doo, P.rs, P.rc)

      if iszero(P.S)
        continue
      end

      e  = 0.0

      # O-O Interactions
      e += _vdw!(F, u, P.Fbuf, P.rbuf, lat, o1, o2, P.ϵ, P.σ; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, o1, o2, P.Qo, P.Qo; S=P.S[1])

      # H-H Coulomb
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, h1, h3, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, h1, h4, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, h2, h3, P.Qh, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, h2, h4, P.Qh, P.Qh; S=P.S[1])

      # H-O Coulomb
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, o1, h3, P.Qo, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, o1, h4, P.Qo, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, o2, h1, P.Qo, P.Qh; S=P.S[1])
      e += _Coulomb!(F, u, P.Fbuf, P.rbuf, lat, o2, h2, P.Qo, P.Qh; S=P.S[1])

      E += P.S[1] * e

      if P.S[2] != 0.0
        P.rbuf2 ./= doo
        P.Fbuf   .= -P.S[2] * e .* P.rbuf2
        F[o1]   .-= P.Fbuf
        F[o2]   .+= P.Fbuf
      end
    end
  end

  E
end