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

    roo = norm(u[o2] - u[o1])
    srhat = (u[o2] - u[o1]) ./ roo
    switchSR!(S, roo, P.rs, P.rc)

    if iszero(S)
      continue
    end

    e = 0.0

    e += _vdw!(F, u, o1, o2, P.ϵoo, P.σoo; S=S[1])

    for i in [h1, h2]
      for j in [h3, h4]
        e += _Coulomb!(F, u, i, j, P.Qh, P.Qh; S=S[1])
      end
    end

    e += _getMforces!(F, u, par[1], par[2], P.drel, P.Qh, P.Qm, S)

    E += S[1] * e

    if S[2] != 0.0
      f .= -S[2] * e .* srhat
      spreadForce!(F, -f, par[1], wa)
      spreadForce!(F,  f, par[2], wa)
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
        E += tip4pf_pbc!(F, u, a, b, S, NC, lat, P)

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
  r1o = u[h1] - u[o1]
  r2o = u[h2] - u[o1]
  r12 = u[h2] - u[h1]
  r3o = u[h3] - u[o2]
  r4o = u[h4] - u[o2]
  r34 = u[h4] - u[h3]

  # Get a
  a1 = 1 / (1 + (norm(r2o) / norm(r1o)))
  a2 = 1 / (1 + (norm(r4o) / norm(r3o)))

  # Get M-sites
  x1 = r1o .+ (a1 * r12)
  m1 = u[o1] + drel * (x1 / norm(x1))
  x2 = r3o .+ (a2 * r34)
  m2 = u[o2] + drel * (x2 / norm(x2))

  # Get Gamma
  γ1 = drel / norm(x1)
  γ2 = drel / norm(x2)

  (a1, a2, γ1, γ2), (m1, m2)
end

function spreadMforces!(
  F::AbstractVector, Fd::AbstractVector, rid::AbstractVector,
  w::Vector{Int64}, γ::Float64, a::Float64
)

  F1 = dot(rid, Fd) / dot(rid, rid) * rid
  X  = Fd - F1

  F[w[1]] .+= Fd - γ * X
  F[w[2]] .+= (1 - a) * γ * X
  F[w[3]] .+= a * γ * X

end

function _getMforces!(
  F::Vector{Af}, u::Vector{Au}, w1::V, w2::V, 
  drel::Fl, Qh::Fl, Qm::Fl, S::AbstractVector
) where {Af <: AbstractVector, Au <: AbstractVector, V <: Vector{Int64}, Fl <: Float64}
  o1, h1, h2 = w1
  o2, h3, h4 = w2

  (a1, a2, γ1, γ2), (m1, m2) = getMsiteVars(u, w1, w2, drel)

  rid1 = m1 - u[o1]
  rid2 = m2 - u[o2]

  # H1 -- M2
  E,f     = _Coulomb(u[h1], m2, Qh, Qm)
  f     .*= S[1]
  F[h1] .-= f
  spreadMforces!(F, f, rid2, w2, γ2, a2)

  # H2 -- M2
  e,f     = _Coulomb(u[h2], m2, Qh, Qm)
  f     .*= S[1]
  E      += e
  F[h2] .-= f
  spreadMforces!(F, f, rid2, w2, γ2, a2)

  # H3 -- M1
  e,f     = _Coulomb(u[h3], m1, Qh, Qm)
  f     .*= S[1]
  E      += e
  F[h3] .-= f
  spreadMforces!(F, f, rid1, w1, γ1, a1)

  # H4 -- M1
  e,f     = _Coulomb(u[h4], m1, Qh, Qm)
  f     .*= S[1]
  E      += e
  F[h4] .-= f
  spreadMforces!(F, f, rid1, w1, γ1, a1)

  # M1 -- M2
  e,f     = _Coulomb(m1, m2, Qm, Qm)
  f     .*= S[1]
  E      += e
  spreadMforces!(F, -f, rid1, w1, γ1, a1)
  spreadMforces!(F,  f, rid2, w2, γ2, a2)

  E
end

function tip4pf_pbc!(
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
        F[o1] .-= f * S[1]
        F[o2] .+= f * S[1]

        # H1-H3 Coulomb
        e,f    = _Coulomb(u[h1], h3t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h1] .-= f * S[1]
        F[h3] .+= f * S[1]

        # H2-H4 Coulomb
        e,f    = _Coulomb(u[h1], h4t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h1] .-= f * S[1]
        F[h4] .+= f * S[1]

        # H2-H3 Coulomb
        e,f    = _Coulomb(u[h2], h3t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h2] .-= f * S[1]
        F[h3] .+= f * S[1]

        # H2-H4 Coulomb
        e,f    = _Coulomb(u[h2], h4t, p.Qh, p.Qh)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        F[h2] .-= f * S[1]
        F[h4] .+= f * S[1]

        # Spread dS force
        spreadForce!(F, -Fs, w1, wa)
        spreadForce!(F,  Fs, w2, wa)
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

  (a1, a2, γ1, γ2), (m1, m2) = getMsiteVars(u, w1, w2, drel)

  rid1 = m1 - u[o1]
  rid2 = m2 - u[o2]

  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        ot  = u[o2] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        roo = norm(ot - u[o1])
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
        f     .*= S[1]
        F[h1] .-= f
        spreadMforces!(F, f, rid2, w2, γ2, a2)

        # H2 -- M2
        e,f     = _Coulomb(u[h2], m2t, Qh, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        f     .*= S[1]
        F[h2] .-= f
        spreadMforces!(F, f, rid2, w2, γ2, a2)

        # H3 -- M1
        e,f     = _Coulomb(h3t, m1, Qh, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        f     .*= S[1]
        F[h3] .-= f
        spreadMforces!(F, f, rid1, w1, γ1, a1)

        # H4 -- M1
        e,f     = _Coulomb(h4t, m1, Qh, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        f     .*= S[1]
        F[h4] .-= f
        spreadMforces!(F, f, rid1, w1, γ1, a1)

        # M1 -- M2
        e,f     = _Coulomb(m1, m2t, Qm, Qm)
        E      += e * S[1]
        Fs    .+= -S[2] * e .* srhat
        f     .*= S[1]
        spreadMforces!(F, -f, rid1, w1, γ1, a1)
        spreadMforces!(F,  f, rid2, w2, γ2, a2)

        # Spread dS force
        spreadForce!(F, -Fs, w1, wa)
        spreadForce!(F,  Fs, w2, wa)
      end
    end
  end

  E
end