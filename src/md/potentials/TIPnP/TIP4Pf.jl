"""
TIP4P/2005f 

Based on:
González, M. A., & Abascal, J. L. (2011). A flexible model for water based on
TIP4P/2005. The Journal of chemical physics, 135(22).

Link:
https://pubs.aip.org/aip/jcp/article/135/22/224516/190786
"""
TIP4Pf(; constraints=nothing) = Calculator(TIP4Pf; EF=TIP4Pf!, constraints=constraints)

struct _TIP4Pf_PotVars{
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
TIP4Pf(bdys::Union{MyCell, Vector{MyAtoms}}) = _TIP4Pf_PotVars(
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