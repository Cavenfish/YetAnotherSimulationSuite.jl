"""
CO-CO Potential

Based on:
van Hemert, Marc C., Junko Takahashi, and Ewine F. van Dishoeck. "Molecular
dynamics study of the photodesorption of CO ice." The Journal of Physical
Chemistry A 119.24 (2015): 6354-6369.

Link:
https://pubs.acs.org/doi/full/10.1021/acs.jpca.5b02611
"""
MvHff(; constraints=nothing) = Calculator(MvHff; EF=MvHff!, constraints=constraints)

struct _MvHff_PotVars{
  F<:Float64, AV2D<:AbstractVector, AV3D<:AbstractVector
} <: PotVars
  D::F
  a::F
  req::F
  Acc::F
  Aco::F
  Aoo::F
  Bcc::F
  Bco::F
  Boo::F
  Ccc::F
  Cco::F
  Coo::F
  Qc::F
  Qo::F
  αc::F
  αo::F
  rs::F
  rc::F
  S::AV2D
  Fbuf::AV3D
  rbuf::AV3D
  r1::AV3D
  r2::AV3D
  x1::AV3D
  x2::AV3D
  F_Q1::AV3D
  F_Q2::AV3D
end

MvHff(bdys::Union{Vector{MyAtoms}, MyCell}) = _MvHff_PotVars(
  11.230139012256362,
  2.3281,
  1.1282058129221093,
  361.367206403597,
  1516.76265699823,
  6370.185468304371,
  2.8345891887553925,
  3.5432364859442407,
  4.2518837831330885,
  33.44955570065988,
  15.189133724736946,
  10.546349734130885,
  -1.7835026375774934,
  -2.333732174702465,
  3.843702939952312,
  2.131611069944055,
  11.0,
  12.0,
  MVector{2}(zeros(2)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3))
)

function MvHff!(F, u, p)
  E   = 0.0
  P   = p.potVars
  lat = isa(p, optVars) ? p.lattice : p.ensemble.lattice

  for i = 1:length(p.mols)
    c1, o1 = p.mols[i]

    E += _Morse!(F, u, lat, c1, o1, P.D, P.a, P.req, P.rc)

    for j = i+1:length(p.mols)
      c2, o2 = p.mols[j]

      doo = pbcVec!(P.rbuf, u[o1], u[o2], P.rc, lat)
      switchLR!(P.S, doo, P.rs, P.rc)

      if iszero(P.S)
        continue
      end

      e = 0.0

      #C--C
      e += _Buckingham!(F, u, lat, c1, c2, P.Acc, P.Bcc, P.Ccc, P.rc; S=P.S[1])

      #C--O
      e += _Buckingham!(F, u, lat, c1, o2, P.Aco, P.Bco, P.Cco, P.rc; S=P.S[1])

      #O--C
      e += _Buckingham!(F, u, lat, o1, c2, P.Aco, P.Bco, P.Cco, P.rc; S=P.S[1])

      #O--O
      e += _Buckingham!(F, u, lat, o1, o2, P.Aoo, P.Boo, P.Coo, P.rc; S=P.S[1])

      #Special Electrostatics
      e += _electroMvH!(F, u, c1, o1, c2, o2, lat, P)

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

function _electroMvH!(F, u, c1, o1, c2, o2, lat, P)

  #Bond lengths
  d1 = pbcVec!(P.r1, u[c1], u[o1], P.rc, lat)
  d2 = pbcVec!(P.r2, u[c2], u[o2], P.rc, lat)

  #CM weights
  wo = 0.57135
  wc = 0.42865

  #X locations
  @. P.x1 = (wc * u[c1]) + (wo * u[o1])
  @. P.x2 = (wc * u[c2]) + (wo * u[o2])

  #Variable charges
  Qc1 = P.Qc * exp(-P.αc * (d1 - P.req))
  Qo1 = P.Qo * exp(-P.αo * (d1 - P.req))
  Qx1 = -(Qc1 + Qo1)
  Qc2 = P.Qc * exp(-P.αc * (d2 - P.req))
  Qo2 = P.Qo * exp(-P.αo * (d2 - P.req))
  Qx2 = -(Qc2 + Qo2)

  #The additional contribution to the forces resulting from
  #the position dependent charges is: -α * Q * rvec / r

  #C--C
  ϵ      = _Coulomb!(P.Fbuf, u[c1], u[c2], lat, Qc1, Qc2, P.rc)
  E      = ϵ
  @. P.F_Q1   = P.αc * ϵ * P.r1 / d1
  @. P.F_Q2   = P.αc * ϵ * P.r2 / d2
  @. F[c1] -= (P.Fbuf + P.F_Q1) * P.S[1]
  @. F[o1] += P.F_Q1 * P.S[1]
  @. F[c2] += P.Fbuf - P.F_Q2 * P.S[1]
  @. F[o2] += P.F_Q2 * P.S[1]

# C-X
  ϵ    = _Coulomb!(P.Fbuf, u[c1], P.x2, lat, Qc1, Qx2, P.rc)
  E     += ϵ
  @. P.F_Q1   = P.αc * ϵ * P.r1 / d1
  @. P.F_Q2   = - (P.αc * Qc2 + P.αo * Qo2) * ϵ/Qx2 * P.r2 / d2
  @. F[c1] -= (P.Fbuf + P.F_Q1) * P.S[1]
  @. F[o1] += P.F_Q1 * P.S[1]
  @. F[c2] += (wc * P.Fbuf) - P.F_Q2 * P.S[1]
  @. F[o2] += (wo * P.Fbuf) + P.F_Q2 * P.S[1]

  # C-O
  ϵ   = _Coulomb!(P.Fbuf, u[c1], u[o2], lat, Qc1, Qo2, P.rc)
  E     += ϵ
  @. P.F_Q1   = P.αc * ϵ * P.r1 / d1
  @. P.F_Q2   = P.αo * ϵ * P.r2 / d2
  @. F[c1] -= (P.Fbuf + P.F_Q1) * P.S[1]
  @. F[o1] += P.F_Q1 * P.S[1]
  @. F[c2] -= P.F_Q2 * P.S[1]
  @. F[o2] += (P.Fbuf + P.F_Q2) * P.S[1]

  # X-C
  ϵ   = _Coulomb!(P.Fbuf, P.x1, u[c2], lat, Qx1, Qc2, P.rc)
  E     += ϵ
  @. P.F_Q1   = - (P.αc * Qc1 + P.αo * Qo1) * ϵ/Qx1 * P.r1 / d1
  @. P.F_Q2   = P.αc * ϵ * P.r2 / d2
  @. F[c1] -= ((wc * P.Fbuf) + P.F_Q1) * P.S[1]
  @. F[o1] += P.F_Q1 - (wo * P.Fbuf) * P.S[1]
  @. F[c2] += (P.Fbuf - P.F_Q2) * P.S[1]
  @. F[o2] += P.F_Q2 * P.S[1]

  # X-X
  ϵ   = _Coulomb!(P.Fbuf, P.x1, P.x2, lat, Qx1, Qx2, P.rc)
  E     += ϵ
  @. P.F_Q1   = - (P.αc * Qc1 + P.αo * Qo1) * ϵ/Qx1 * P.r1 / d1
  @. P.F_Q2   = - (P.αc * Qc2 + P.αo * Qo2) * ϵ/Qx2 * P.r2 / d2
  @. F[c1] += ((-wc * P.Fbuf) - P.F_Q1) * P.S[1]
  @. F[o1] += ((-wo * P.Fbuf) + P.F_Q1) * P.S[1]
  @. F[c2] += (( wc * P.Fbuf) - P.F_Q2) * P.S[1]
  @. F[o2] += (( wo * P.Fbuf) + P.F_Q2) * P.S[1]


  # X-O
  ϵ   = _Coulomb!(P.Fbuf, P.x1, u[o2], lat, Qx1, Qo2, P.rc)
  E     += ϵ
  @. P.F_Q1   = - (P.αc * Qc1 + P.αo * Qo1) * ϵ/Qx1 * P.r1 / d1
  @. P.F_Q2   = P.αo * ϵ * P.r2 / d2
  @. F[c1] -= ((wc * P.Fbuf) + P.F_Q1) * P.S[1]
  @. F[o1] += (P.F_Q1 - (wo * P.Fbuf)) * P.S[1]
  @. F[c2] -= P.F_Q2 * P.S[1]
  @. F[o2] += (P.Fbuf + P.F_Q2) * P.S[1]

  # O-C
  ϵ   = _Coulomb!(P.Fbuf, u[o1], u[c2], lat, Qo1, Qc2, P.rc)
  E     += ϵ
  @. P.F_Q1   = P.αo * ϵ * P.r1 / d1
  @. P.F_Q2   = P.αc * ϵ * P.r2 / d2
  @. F[c1] -= P.F_Q1 * P.S[1]
  @. F[o1] += (P.F_Q1 - P.Fbuf) * P.S[1]
  @. F[c2] += (P.Fbuf - P.F_Q2) * P.S[1]
  @. F[o2] += P.F_Q2 * P.S[1]

  # O-X
  ϵ   = _Coulomb!(P.Fbuf, u[o1], P.x2, lat, Qo1, Qx2, P.rc)
  E     += ϵ
  @. P.F_Q1   = P.αo * ϵ * P.r1 / d1
  @. P.F_Q2   = - (P.αc * Qc2 + P.αo * Qo2) * ϵ/Qx2 * P.r2 / d2
  @. F[c1] -= P.F_Q1 * P.S[1]
  @. F[o1] += (P.F_Q1 - P.Fbuf) * P.S[1]
  @. F[c2] += ((wc * P.Fbuf) - P.F_Q2) * P.S[1]
  @. F[o2] += ((wo * P.Fbuf) + P.F_Q2) * P.S[1]

  # O-O
  ϵ   = _Coulomb!(P.Fbuf, u[o1], u[o2], lat, Qo1, Qo2, P.rc)
  E     += ϵ
  @. P.F_Q1   = P.αo * ϵ * P.r1 / d1
  @. P.F_Q2   = P.αo * ϵ * P.r2 / d2
  @. F[c1] -= P.F_Q1 * P.S[1]
  @. F[o1] += (P.F_Q1 - P.Fbuf) * P.S[1]
  @. F[c2] -= P.F_Q2 * P.S[1]
  @. F[o2] += (P.Fbuf + P.F_Q2) * P.S[1]
  

  E
end