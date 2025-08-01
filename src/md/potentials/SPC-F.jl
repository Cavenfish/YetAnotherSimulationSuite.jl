"""
SPC-F
"""

SPCF(; constraints=nothing) = Calculator(SPCF; EF=SPCF!, constraints=constraints)

struct _SPCF_PotVars{F<:Float64} <: PotVars
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
  10.0
)

function SPCF!(F, u, p)
  E  = 0.0
  wa = [-1.0, 0.5, 0.5]
  f  = zeros(3)
  S  = zeros(2)
  P  = p.potVars
  NC = p.PBC .* p.NC

  for mol in p.mols
    o, h1, h2 = mol

    E += _harmonicBond!(F, u, o, h1, P.Kb, P.req)
    E += _harmonicBond!(F, u, o, h2, P.Kb, P.req)
    E += _harmonicBondAngle!(F, u, h1, o, h2, P.Kθ, P.θeq)
  end

  for par in p.pars
    o1, h1, h2 = par[1]
    o2, h3, h4 = par[2]

    roo = norm(u[o1] - u[o2])
    srhat = (u[o2] - u[o1]) / roo
    switchLR!(S, roo, P.rs, P.rc)

    if iszero(S)
      continue
    end

    e  = 0.0

    e += _vdw!(F, u, o1, o2, P.ϵ, P.σ)
    e += _Coulomb!(F, u, o1, o2, P.Qo, P.Qo)

    for i in [h1, h2]
      for j in [h3, h4]
        e += _Coulomb!(F, u, i, j, P.Qh, P.Qh)
      end
    end

    for i in [h1,h2]
      e += _Coulomb!(F, u, o2, i, P.Qo, P.Qh)
    end

    for i in [h3,h4]
      e += _Coulomb!(F, u, o1, i, P.Qo, P.Qh)
    end

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
      
      E += spcf_selfInter_pbc!(F, u, a, S, NC, lat, P)

      for j = i+1:length(p.mols)
        b = p.mols[j]
      
        E += spcf_inter_pbc!(F, u, a, b, S, NC, lat, P)

      end
    end

  end

  E
end

function spcf_selfInter_pbc!(
  F::Vector{Af}, u::Vector{Au}, w::Vector{Int64}, S::AbstractVector,
  NC::Vector{Int64}, L::AbstractMatrix, p::P
) where {Af, Au, P}

  E = 0.0
  wa = [-1.0, 0.5, 0.5]
  Fs = zeros(3)
  o, h1, h2 = w

  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        ot  = u[o] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        roo = norm(u[o] - ot)
        srhat = (ot - u[o]) / roo
        switchLR!(S, roo, p.rs, p.rc)

        if iszero(S)
          continue
        end

        h1t = u[h1] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        h2t = u[h2] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        if i >= 0 && j >= 0 && k >= 0
          # O-O VDW
          e,f    = _vdw(u[o], ot, p.ϵ, p.σ)
          E     += e * S[1]
          Fs    .= -S[2] * e .* srhat
          F[o] .-= f
          F[o] .+= f

          # O-O Coulomb
          e,f    = _Coulomb(u[o], ot, p.Qo, p.Qo)
          E      += e * S[1]
          Fs   .+= -S[2] * e .* srhat
          F[o] .-= f
          F[o] .+= f

          # H1-H1 Coulomb
          e,f    = _Coulomb(u[h1], h1t, p.Qh, p.Qh)
          E      += e * S[1]
          Fs   .+= -S[2] * e .* srhat
          F[h1] .-= f
          F[h1] .+= f

          # H2-H2 Coulomb
          e,f    = _Coulomb(u[h2], h2t, p.Qh, p.Qh)
          E      += e * S[1]
          Fs   .+= -S[2] * e .* srhat
          F[h2] .-= f
          F[h2] .+= f
        end

        # O-H1 Coulomb
        e,f    = _Coulomb(u[o], h1t, p.Qo, p.Qh)
        E      += e * S[1]
        Fs   .+= -S[2] * e .* srhat
        F[o]  .-= f
        F[h1] .+= f

        # O-H2 Coulomb
        e,f    = _Coulomb(u[o], h2t, p.Qo, p.Qh)
        E      += e * S[1]
        Fs   .+= -S[2] * e .* srhat
        F[o]  .-= f
        F[h2] .+= f

        # Spread dS force
        spreadForce!(F,  Fs, w, wa)
        spreadForce!(F, -Fs, w, wa) 
      end
    end
  end

  E
end

function spcf_inter_pbc!(
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

        o2t = u[o2] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        roo = norm(u[o1] - o2t)
        srhat = (o2t - u[o1]) / roo
        switchLR!(S, roo, p.rs, p.rc)

        if iszero(S)
          continue
        end

        h3t = u[h3] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)
        h4t = u[h4] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        # O-O VDW
        e,f    = _vdw(u[o1], o2t, p.ϵ, p.σ)
        E     += e * S[1]
        Fs     .= -S[2] * e * srhat
        F[o1] .-= f
        F[o2] .+= f

        # O-O Coulomb
        e,f    = _Coulomb(u[o1], o2t, p.Qo, p.Qo)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[o1] .-= f
        F[o2] .+= f

        # H1-H3 Coulomb
        e,f    = _Coulomb(u[h1], h3t, p.Qh, p.Qh)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[h1] .-= f
        F[h3] .+= f

        # H1-H4 Coulomb
        e,f    = _Coulomb(u[h1], h4t, p.Qh, p.Qh)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[h1] .-= f
        F[h4] .+= f

        # H2-H3 Coulomb
        e,f    = _Coulomb(u[h2], h3t, p.Qh, p.Qh)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[h2] .-= f
        F[h3] .+= f

        # H2-H4 Coulomb
        e,f    = _Coulomb(u[h2], h4t, p.Qh, p.Qh)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[h2] .-= f
        F[h4] .+= f

        # H1-O2 Coulomb
        e,f    = _Coulomb(u[h1], o2t, p.Qh, p.Qo)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[h1] .-= f
        F[o2] .+= f

        # H2-O2 Coulomb
        e,f    = _Coulomb(u[h2], o2t, p.Qh, p.Qo)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[h2] .-= f
        F[o2] .+= f

        # H3-O1 Coulomb
        e,f    = _Coulomb(u[o1], h3t, p.Qo, p.Qh)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[o1] .-= f
        F[h3] .+= f

        # H4-O1 Coulomb
        e,f    = _Coulomb(u[o1], h4t, p.Qo, p.Qh)
        E      += e
        Fs    .+= -S[2] * e * srhat
        F[o1] .-= f
        F[h4] .+= f

        # Spread dS force
        spreadForce!(F,  Fs, w1, wa)
        spreadForce!(F, -Fs, w2, wa) 
      end
    end
  end

  E
end