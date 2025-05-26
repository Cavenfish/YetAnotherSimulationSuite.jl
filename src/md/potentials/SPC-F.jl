"""
SPC-F
"""

SPCF() = Calculator(SPCF; EF=SPCF!)

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

SPCF(bdys::Vector{MyAtoms}) = _SPCF_PotVars(
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
  E = 0.0
  P = p.potVars

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

  E
end