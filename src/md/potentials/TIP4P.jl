"""
TIP4P/2005f 
"""

TIP4Pf(; constraints=nothing) = Calculator(TIP4Pf; EF=TIP4Pf!, constraints=constraints)

struct _TIP4P_PotVars{F<:Float64} <: PotVars
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
  -2 * 2.1113635
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

  for par in p.pars
    o1, h1, h2 = par[1]
    o2, h3, h4 = par[2]

    E += _vdw!(F, u, o1, o2, P.ϵoo, P.σoo)

    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(F, u, i, j, P.Qh, P.Qh)
      end
    end

    E += _getMforces!(F, u, par[1], par[2], P.drel, P.Qh, P.Qm)
  end

  if any(p.PBC)
    NC    = p.NC .* p.PBC
    all_o = collect(1:3:length(u))
    all_h = [i for i = 1:length(u) if !(i in all_o)]

    for i = 1:length(all_o)
      a = all_o[i]

      E += _pbc!(F, u, a, a, _vdw, p.lattice, NC, (P.ϵoo, P.σoo); cutoff=45.0) / 2
      for j = i+1:length(all_o)
        b = all_o[j]

        E += _pbc!(F, u, a, b, _vdw, p.lattice, NC, (P.ϵoo, P.σoo); cutoff=45.0)
      end
    end

    for i = 1:length(all_h)
      a = all_h[i]

      E += _pbc!(F, u, a, a, _Coulomb, p.lattice, NC, (P.Qh, P.Qh); cutoff=45.0) / 2
      for j = i+1:length(all_h)
        b = all_h[j]

        E += _pbc!(F, u, a, b, _Coulomb, p.lattice, NC, (P.Qh, P.Qh); cutoff=45.0)
      end
    end

    for i = 1:length(p.mols)
      a = p.mols[i]

      E += pbc_Mforces!(F, u, a, a, P.drel, P.Qh, P.Qm, NC, p.lattice; cutoff=45.0) / 2
      for j = i+1:length(p.mols)
        b = p.mols[j]

        E += pbc_Mforces!(F, u, a, b, P.drel, P.Qh, P.Qm, NC, p.lattice; cutoff=45.0)
      end
    end
  end

  E
end

function getMsiteVars(
  u::Vector{A}, w1::V, w2::V, drel::Float64
) where {V <: Vector{Int64}, A <: AbstractVector}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  # Get r vectors
  r1o  = u[h1] - u[o1]
  r2o  = u[h2] - u[o1]
  r3o  = u[h3] - u[o2]
  r4o  = u[h4] - u[o2]

  # Get Angle
  θ1  = getAngle(r1o, r2o)
  θ2  = getAngle(r3o, r4o)

  # Get Norms
  d1o  = norm(r1o)
  d2o  = norm(r2o)
  d3o  = norm(r3o)
  d4o  = norm(r4o)

  # Get dm
  dm1 = drel * (d1o*cos(θ1/2) + d2o*cos(θ1/2))
  dm2 = drel * (d3o*cos(θ2/2) + d4o*cos(θ2/2))

  # Get bisectors
  rbi1 = r1o/d1o + r2o/d2o
  rbi2 = r3o/d3o + r4o/d4o

  # Get weights
  wh1 = dm1 / d1o / norm(rbi1)
  wh2 = dm1 / d2o / norm(rbi1)
  wh3 = dm2 / d3o / norm(rbi2)
  wh4 = dm2 / d4o / norm(rbi2)

  # Get m site vectors
  m1 = u[o1] + wh1*r1o + wh2*r2o
  m2 = u[o2] + wh3*r3o + wh4*r4o

  (wh1, wh2, wh3, wh4), (m1, m2)
end

function _getMforces!(
  F::Vector{Af}, u::Vector{Au}, w1::V, w2::V, drel::Fl, Qh::Fl, Qm::Fl
) where {Af <: AbstractVector, Au <: AbstractVector, V <: Vector{Int64}, Fl <: Float64}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  (wh1, wh2, wh3, wh4), (m1, m2) = getMsiteVars(u, w1, w2, drel)

  # H1 -- M2
  E,f     = _Coulomb(u[h1], m2, Qh, Qm)
  F[h1] .-= f
  F[h3] .+= f * wh3
  F[h4] .+= f * wh4
  F[o2] .+= f * (1 - wh3 - wh4)

  # H2 -- M2
  e,f     = _Coulomb(u[h2], m2, Qh, Qm)
  E      += e
  F[h2] .-= f
  F[h3] .+= f * wh3
  F[h4] .+= f * wh4
  F[o2] .+= f * (1 - wh3 - wh4)

  # H3 -- M1
  e,f     = _Coulomb(u[h3], m1, Qh, Qm)
  E      += e
  F[h3] .-= f
  F[h1] .+= f * wh1
  F[h2] .+= f * wh2
  F[o1] .+= f * (1 - wh1 - wh2)

  # H4 -- M1
  e,f     = _Coulomb(u[h4], m1, Qh, Qm)
  E      += e
  F[h4] .-= f
  F[h1] .+= f * wh1
  F[h2] .+= f * wh2
  F[o1] .+= f * (1 - wh1 - wh2)

  # M1 -- M2
  e,f     = _Coulomb(m1, m2, Qm, Qm)
  E      += e
  F[h1] .-= f * wh1
  F[h2] .-= f * wh2
  F[o1] .-= f * (1 - wh1 - wh2)
  F[h3] .+= f * wh3
  F[h4] .+= f * wh4
  F[o2] .+= f * (1 - wh3 - wh4)

  E
end

function pbc_Mforces!(
  F::Vector{Af}, u::Vector{Au}, w1::V, w2::V, 
  drel::Fl, Qh::Fl, Qm::Fl, NC::V, L::AbstractMatrix;
  cutoff=20.0
) where {Af, Au, V <: Vector{Int64}, Fl <: Float64}
  E = 0.0
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  # buffers
  t   = zeros(3)
  m2t = zeros(3)
  h3t = zeros(3)
  h4t = zeros(3)

  (wh1, wh2, wh3, wh4), (m1, m2) = getMsiteVars(u, w1, w2, drel)

  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        t   .= (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        m2t .= m2 + t
        h3t .= u[h3] + t
        h4t .= u[h4] + t

        norm(m1 - m2t) > cutoff && continue

        # H1 -- M2
        e,f     = _Coulomb(u[h1], m2t, Qh, Qm)
        E      += e
        F[h1] .-= f
        F[h3] .+= f * wh3
        F[h4] .+= f * wh4
        F[o2] .+= f * (1 - wh3 - wh4)

        # H2 -- M2
        e,f     = _Coulomb(u[h2], m2t, Qh, Qm)
        E      += e
        F[h2] .-= f
        F[h3] .+= f * wh3
        F[h4] .+= f * wh4
        F[o2] .+= f * (1 - wh3 - wh4)

        # H3 -- M1
        e,f     = _Coulomb(h3t, m1, Qh, Qm)
        E      += e
        F[h3] .-= f
        F[h1] .+= f * wh1
        F[h2] .+= f * wh2
        F[o1] .+= f * (1 - wh1 - wh2)

        # H4 -- M1
        e,f     = _Coulomb(h4t, m1, Qh, Qm)
        E      += e
        F[h4] .-= f
        F[h1] .+= f * wh1
        F[h2] .+= f * wh2
        F[o1] .+= f * (1 - wh1 - wh2)

        # M1 -- M2
        e,f     = _Coulomb(m1, m2t, Qm, Qm)
        E      += e
        F[h1] .-= f * wh1
        F[h2] .-= f * wh2
        F[o1] .-= f * (1 - wh1 - wh2)
        F[h3] .+= f * wh3
        F[h4] .+= f * wh4
        F[o2] .+= f * (1 - wh3 - wh4)
      end
    end
  end

  E
end