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
  2.287,      # \AA^-1
  0.9419,     # \AA
  3.81209321, # eV rad^-2  1.16123e-3 * (360/2pi)^2
  1.87448,    # rad        107.4  * (2pi/360)
  8.03e-3,    # eV
  3.1644,     # \AA
  0.1546,    # \AA
  2.1113635,  #  qh * Hartree^0.5 * Bohr^0.5
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
  E   = 0.0
  P   = p.potVars
  lat = isa(p, optVars) ? p.lattice : p.ensemble.lattice

  for i = 1:length(p.mols)
    o1, h1, h2 = p.mols[i]

    E += morse!(F, u, lat, o1, h1, P.D, P.a, P.req, P.rc)
    E += morse!(F, u, lat, o1, h2, P.D, P.a, P.req, P.rc)
    E += harmonicBondAngle!(F, u, lat, h1, o1, h2, P.K, P.θeq, P.rc)

    for j = i+1:length(p.mols)
      o2, h3, h4 = p.mols[j]

      doo = pbcVec!(P.rbuf, u[o1], u[o2], P.rc, lat)
      switchLR!(P.S, doo, P.rs, P.rc)

      if iszero(P.S)
        continue
      end

      e = 0.0

      # O-O vdw Interaction
      e += vdw!(F, u, lat, o1, o2, P.ϵoo, P.σoo, P.rc; S=P.S[1])

      # H-H Coulomb Interactions
      e += coulomb!(F, u, lat, h1, h3, P.Qh, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, h1, h4, P.Qh, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, h2, h3, P.Qh, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, h2, h4, P.Qh, P.Qh, P.rc; S=P.S[1])

      e += getMforces!(F, u, p.mols[i], p.mols[j], lat, P)

      E += P.S[1] * e

      if P.S[2] != 0.0
        P.rbuf ./= doo
        P.Fbuf  .= -P.S[2] * e .* P.rbuf
        F[o1]  .-= P.Fbuf
        F[o2]  .+= P.Fbuf
      end
    end
  end

  E
end

function getMsiteVars!(
  P::T, u::Vu, w1::V, w2::V, lat::AbstractMatrix
) where {T, V<:Vector{Int64}, Vu<:AbstractVector}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  # Get r vectors
  _ = pbcVec!(P.r1o, u[o1], u[h1], P.rc, lat)
  _ = pbcVec!(P.r2o, u[o1], u[h2], P.rc, lat)
  _ = pbcVec!(P.r12, u[h1], u[h2], P.rc, lat)
  _ = pbcVec!(P.r3o, u[o2], u[h3], P.rc, lat)
  _ = pbcVec!(P.r4o, u[o2], u[h4], P.rc, lat)
  _ = pbcVec!(P.r34, u[h3], u[h4], P.rc, lat)

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

function getMforces!(
  F::Vf, u::Vu, w1::V, w2::V, lat::AbstractMatrix, P::T
) where {Vf<:AbstractVector, Vu<:AbstractVector, V<:Vector{Int64}, T}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  γ1, γ2 = getMsiteVars!(P, u, w1, w2, lat)

  _ = pbcVec!(P.r1o, u[o1], P.m1, P.rc, lat)
  _ = pbcVec!(P.r2o, u[o2], P.m2, P.rc, lat)

  # H1 -- M2
  E        = coulomb!(P.Fbuf, u[h1], P.m2, lat, P.Qh, P.Qm, P.rc)
  P.Fbuf .*= P.S[1]
  F[h1]  .-= P.Fbuf
  spreadMforces!(F, P.Fbuf, P.r2o, P, w2, γ2)

  # H2 -- M2
  e        = coulomb!(P.Fbuf, u[h2], P.m2, lat, P.Qh, P.Qm, P.rc)
  P.Fbuf .*= P.S[1]
  E       += e
  F[h2]  .-= P.Fbuf
  spreadMforces!(F, P.Fbuf, P.r2o, P, w2, γ2)

  # H3 -- M1
  e        = coulomb!(P.Fbuf, u[h3], P.m1, lat, P.Qh, P.Qm, P.rc)
  P.Fbuf .*= P.S[1]
  E       += e
  F[h3]  .-= P.Fbuf
  spreadMforces!(F, P.Fbuf, P.r1o, P, w1, γ1)

  # H4 -- M1
  e        = coulomb!(P.Fbuf, u[h4], P.m1, lat, P.Qh, P.Qm, P.rc)
  P.Fbuf .*= P.S[1]
  E       += e
  F[h4]  .-= P.Fbuf
  spreadMforces!(F, P.Fbuf, P.r1o, P, w1, γ1)

  # M1 -- M2
  e        = coulomb!(P.Fbuf, P.m1, P.m2, lat, P.Qm, P.Qm, P.rc)
  P.Fbuf .*= P.S[1]
  E       += e
  spreadMforces!(F, P.Fbuf, P.r2o, P, w2, γ2)
  P.Fbuf .*= -1
  spreadMforces!(F, P.Fbuf, P.r1o, P, w1, γ1)

  E
end