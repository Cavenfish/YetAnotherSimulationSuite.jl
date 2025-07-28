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
end

SPCF(bdys::Union{Vector{MyAtoms}, MyCell}) = _SPCF_PotVars(
  48.05913,
  1.0,
  3.97,
  1.910611,
  3.145, 
  0.007,
  -2.959855,
  0.5 * 2.959855
)

function SPCF!(F, u, p)
  E  = 0.0
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

    E += _vdw!(F, u, o1, o2, P.ϵ, P.σ)
    E += _Coulomb!(F, u, o1, o2, P.Qo, P.Qo)


    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(F, u, i, j, P.Qh, P.Qh)
      end
    end

    for i in [h1,h2]
      E += _Coulomb!(F, u, o2, i, P.Qo, P.Qh)
    end

    for i in [h3,h4]
      E += _Coulomb!(F, u, o1, i, P.Qo, P.Qh)
    end
  end

  if any(p.PBC)
    NC    = p.NC .* p.PBC
    lat   = isa(p, optVars) ? p.lattice : p.ensemble.lattice
    all_o = collect(1:3:length(u))
    all_h = [i for i = 1:length(u) if !(i in all_o)]

    for i = 1:length(all_o)
      for j = i+1:length(all_o)
        a = all_o[i]
        b = all_o[j]

        E += _pbc!(F, u, a, b, _vdw, p.lattice, NC, (P.ϵ, P.σ); cutoff=45.0)
        E += _pbc!(F, u, a, b, _Coulomb, p.lattice, NC, (P.Qo, P.Qh); cutoff=45.0)
      end
    end

    for i = 1:length(all_h)
      for j = i+1:length(all_h)
        a = all_h[i]
        b = all_h[j]

        E += _pbc!(F, u, a, b, _Coulomb, p.lattice, NC, (P.Qh, P.Qh); cutoff=45.0)
      end
    end

    for i = 1:length(all_h)
      for j = 1:length(all_o)
        a = all_h[i]
        b = all_o[j]

        E += _pbc!(F, u, a, b, _Coulomb, p.lattice, NC, (P.Qh, P.Qo); cutoff=45.0)
      end
    end

  end

  E
end