"""
TIP4P/2005f 

Based on:
González, M. A., & Abascal, J. L. (2011). A flexible model for water based on
TIP4P/2005. The Journal of chemical physics, 135(22).

Link:
https://pubs.aip.org/aip/jcp/article/135/22/224516/190786
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
  rs::F
  rc::F
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
)

function TIP4Pf!(F, u, p)
  E = 0.0
  P = p.potVars
  S = zeros(2)
  wa = [-1.0, 0.5, 0.5]
  f  = zeros(3)

  for mol in p.mols
    o, h1, h2 = mol

    E += _Morse!(F, u, o, h1, P.D, P.a, P.req)
    E += _Morse!(F, u, o, h2, P.D, P.a, P.req)
    E += _harmonicBondAngle!(F, u, h1, o, h2, P.K, P.θeq)
  end

  for par in p.pars
    o1, h1, h2 = par[1]
    o2, h3, h4 = par[2]

    roo = norm(u[o1] - u[o2])
    srhat = (u[o2] - u[o1]) ./ roo
    switchSR!(S, roo, P.rs, P.rc)

    if iszero(S)
      continue
    end

    e = 0.0

    e += _vdw!(F, u, o1, o2, P.ϵoo, P.σoo)

    for i in [h1, h2]
      for j in [h3, h4]
        e += _Coulomb!(F, u, i, j, P.Qh, P.Qh)
      end
    end

    e += _getMforces!(F, u, par[1], par[2], P.drel, P.Qh, P.Qm)

    E += S[1] * e

    if S[2] != 0.0
      f .= -S[2] * e .* srhat
      spreadForce!(F,  f, par[1], wa)
      spreadForce!(F, -f, par[2], wa)
    end
  end

  if any(p.PBC)
    NC    = p.NC .* p.PBC
    lat   = isa(p, optVars) ? p.lattice : p.ensemble.lattice

    for i = 1:length(p.mols)
      a = p.mols[i]

      # vdw and H-H Coulomb
      E += tip4pf_selfInter_pbc!(F, u, a, S, NC, lat, P)

      # M-site Coulomb
      E += pbc_Mforces!(F, u, a, a, P.drel, P.Qh, P.Qm, NC, lat, P.rs, P.rc, S) / 2

      for j = i+1:length(p.mols)
        b = p.mols[j]

        # vdw and H-H Coulomb
        E += tip4pf_inter_pbc!(F, u, a, b, S, NC, lat, P)

        # M-site Coulomb
        E += pbc_Mforces!(F, u, a, b, P.drel, P.Qh, P.Qm, NC, lat, P.rs, P.rc, S)
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
  F::Vector{Af}, u::Vector{Au}, w1::V, w2::V, 
  drel::Fl, Qh::Fl, Qm::Fl
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

function tip4pf_selfInter_pbc!(
  F::Vector{Af}, u::Vector{Au}, w::Vector{Int64}, S::AbstractVector,
  NC::Vector{Int64}, L::AbstractMatrix, p::P
) where {Af, Au, P}

  E = 0.0
  wa = [-1.0, 0.5, 0.5]
  Fs = zeros(3)
  o, h1, h2 = w

  for i = 0:NC[1]
    for j = 0:NC[2]
      for k = 0:NC[3]
        (i,j,k) == (0,0,0) && continue

        ot  = u[o] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        roo = norm(u[o] - ot)
        srhat = (ot - u[o]) ./ roo
        switchSR!(S, roo, p.rs, p.rc)

        if iszero(S)
          continue
        end

        h1t = u[h1] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        h2t = u[h2] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        # O-O VDW
        e,f    = _vdw(u[o], ot, p.ϵoo, p.σoo)
        E     += e * S[1]
        Fs    .= -S[2] * e .* srhat
        F[o] .-= f
        F[o] .+= f

        # H1-H1 Coulomb
        e,f    = _Coulomb(u[h1], h1t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h1] .-= f
        F[h1] .+= f

        # H2-H2 Coulomb
        e,f    = _Coulomb(u[h2], h2t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h2] .-= f
        F[h2] .+= f

        # H1-H2 Coulomb
        e,f    = _Coulomb(u[h1], h2t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h1] .-= f
        F[h2] .+= f

        # H2-H1 Coulomb (needed since loops are 0:NC not -NC:NC)
        e,f    = _Coulomb(u[h2], h1t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h2] .-= f
        F[h1] .+= f

        # Spread dS force
        spreadForce!(F,  Fs, w, wa)
        spreadForce!(F, -Fs, w, wa)
      end
    end
  end      

  E
end

function tip4pf_inter_pbc!(
  F::Vector{Af}, u::Vector{Au}, w1::Vi, w2::Vi, S::AbstractVector,
  NC::Vector{Int64}, L::AbstractMatrix, p::P
) where {Af, Au, Vi <: Vector{Int64}, P}

  E = 0.0
  wa = [-1.0, 0.5, 0.5]
  Fs = zeros(3)
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        ot  = u[o2] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        roo = norm(u[o1] - ot)
        srhat = (ot - u[o1]) ./ roo
        switchSR!(S, roo, p.rs, p.rc)

        if iszero(S)
          continue
        end

        h3t = u[h3] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        h4t = u[h4] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        # O1-O2 VDW
        e,f    = _vdw(u[o1], ot, p.ϵoo, p.σoo)
        E      += e * S[1]
        Fs     .= -S[2] * e .* srhat
        F[o1] .-= f
        F[o2] .+= f

        # H1-H3 Coulomb
        e,f    = _Coulomb(u[h1], h3t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h1] .-= f
        F[h3] .+= f

        # H2-H4 Coulomb
        e,f    = _Coulomb(u[h1], h4t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h1] .-= f
        F[h4] .+= f

        # H2-H3 Coulomb
        e,f    = _Coulomb(u[h2], h3t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h2] .-= f
        F[h3] .+= f

        # H2-H4 Coulomb
        e,f    = _Coulomb(u[h2], h4t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h2] .-= f
        F[h4] .+= f

        # Spread dS force
        spreadForce!(F,  Fs, w1, wa)
        spreadForce!(F, -Fs, w2, wa)
      end
    end
  end      

  E
end

function pbc_Mforces!(
  F::Vector{Af}, u::Vector{Au}, w1::V, w2::V, 
  drel::Fl, Qh::Fl, Qm::Fl, NC::V, L::AbstractMatrix, 
  rs::Fl, rc::Fl, S::AbstractVector
) where {Af, Au, V <: Vector{Int64}, Fl <: Float64}
  E = 0.0
  wa = [-1.0, 0.5, 0.5]
  Fs = zeros(3)
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

        ot  = u[o2] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        roo = norm(u[o1] - ot)
        srhat = (ot - u[o1]) ./ roo
        switchSR!(S, roo, rs, rc)

        if iszero(S)
          continue
        end

        t   .= (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        m2t .= m2 + t
        h3t .= u[h3] + t
        h4t .= u[h4] + t

        # H1 -- M2
        e,f     = _Coulomb(u[h1], m2t, Qh, Qm)
        E      += e * S[1]
        Fs     .= -S[2] * e .* srhat
        F[h1] .-= f
        F[h3] .+= f * wh3
        F[h4] .+= f * wh4
        F[o2] .+= f * (1 - wh3 - wh4)

        # H2 -- M2
        e,f     = _Coulomb(u[h2], m2t, Qh, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h2] .-= f
        F[h3] .+= f * wh3
        F[h4] .+= f * wh4
        F[o2] .+= f * (1 - wh3 - wh4)

        # H3 -- M1
        e,f     = _Coulomb(h3t, m1, Qh, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h3] .-= f
        F[h1] .+= f * wh1
        F[h2] .+= f * wh2
        F[o1] .+= f * (1 - wh1 - wh2)

        # H4 -- M1
        e,f     = _Coulomb(h4t, m1, Qh, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h4] .-= f
        F[h1] .+= f * wh1
        F[h2] .+= f * wh2
        F[o1] .+= f * (1 - wh1 - wh2)

        # M1 -- M2
        e,f     = _Coulomb(m1, m2t, Qm, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h1] .-= f * wh1
        F[h2] .-= f * wh2
        F[o1] .-= f * (1 - wh1 - wh2)
        F[h3] .+= f * wh3
        F[h4] .+= f * wh4
        F[o2] .+= f * (1 - wh3 - wh4)

        # Spread dS force
        spreadForce!(F,  Fs, w1, wa)
        spreadForce!(F, -Fs, w2, wa)
      end
    end
  end

  E
end