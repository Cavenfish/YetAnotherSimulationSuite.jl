"""
CH4 TTM-nrg Potential
from Riera et al. 2020
https://doi.org/10.1021/acs.jpcb.0c08728
"""

struct _CH4_PotVars{F<:Float64} <: PotVars
  D::F
  a::F
  req::F
  K::F
  θeq::F
  Acc::F
  Ach::F
  Ahh::F
  Bcc::F
  Bch::F
  Bhh::F
  Ccc::F
  Cch::F
  Chh::F
  Qh::F
  Qc::F
  αh::F
  αc::F
  A6cc::F
  A6ch::F
  A6hh::F
  Eq::Vector
  μ::Vector
  α::Vector
  Q::Vector
end

function CH4(bdys::Vector{Atom})
  N  = div(length(bdys), 5)
  Qh = 0.51163
  Qc = -4Qh
  αh = 0.38978
  αc = 1.39327
  μ  = [zeros(3) for i in bdys]
  α  = repeat([αc, αh, αh, αh, αh], N)
  Q  = repeat([Qc, Qh, Qh, Qh, Qh], N)
  Eq = [zeros(3) for i = 1:length(α)]

  _CH4_PotVars(
    4.0017,  # eV
    2.0415,  # \AA^-1
    1.0887,  # \AA
    3.33487, # eV rad^-1   76.9061 * 0.043363
    1.87411, # rad 107.379 * (2pi/360)
    1852.25, #
    141.317, #
    112.503, #
    3.37925, # \AA^-1
    3.25885, # \AA^-1
    4.05972, # \AA^-1
    13.15,   #
    4.51455, #
    1.59497, #
    Qh,      #
    Qc,      #
    αh,      # \AA^3
    αc,      # \AA^3
    αc^(2/6),
    (αc*αh)^(1/6),
    αh^(2/6),
    Eq,
    μ,
    α,
    Q
  )
end

function CH4(dv, v, u, p, t)

  # initialize things
  E = 0.0
  F = zero(u)
  P = p.potVars

  _getPermanentEfield!(P.Eq, u, P.Q, P.α)
  _getDipoles4TTM_MatrixInversion!(P.μ, u, P.Q, P.α, P.Eq)

  for i = 1:length(u)
    E -= 0.5 * dot(P.μ[i], P.Eq[i])
  end

  for mol in p.mols
    c, hs = mol[1], mol[2:5]

    for i in hs
      E += _Morse!(F, u, c, i, P.D, P.a, P.req)
    end

    for i in 1:4
      for j in i+1:4
        E += _harmonicBondAngle!(F, u, hs[i], c, hs[j], P.K, P.θeq)
      end
    end
    
  end

  for par in p.pars
    c1, hs1 = par[1][1], par[1][2:5]
    c2, hs2 = par[2][1], par[2][2:5] 

    E += _interTTM!(
      F, u, P.μ, c1, c2, P.Qc, P.Qc, 
      P.Acc, P.Bcc, P.Ccc, P.A6cc; damp=tangToennies, p=P.Bcc
    )

    k = true
    for i in hs1
      E += _interTTM!(
        F, u, P.μ, i, c2, P.Qh, P.Qc, 
        P.Ach, P.Bch, P.Cch, P.A6ch; damp=tangToennies, p=P.Bch
      )
      for j in hs2
        E += _interTTM!(
          F, u, P.μ, i, j, P.Qh, P.Qh, 
          P.Ahh, P.Bhh, P.Chh, P.A6hh; damp=tangToennies, p=P.Bhh
        )
        if k
          E += _interTTM!(
            F, u, P.μ, j, c1, P.Qh, P.Qc, 
            P.Ach, P.Bch, P.Cch, P.A6ch; damp=tangToennies, p=P.Bch
          )
        end

      end
      k = false
    end
  end
  

  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end

function CH4(F, G, y0, p)

  # initialize things
  P      = p.potVars
  E      = 0.0
  u      = Vector[]
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

  _getPermanentEfield!(P.Eq, u, P.Q, P.α)
  # _getDipoles4TTM_Iterative!(P.μ, u, P.Q, P.α, p.mols, P.Eq)
  _getDipoles4TTM_MatrixInversion!(P.μ, u, P.Q, P.α, P.Eq)

  for i = 1:length(u)
    E -= 0.5 * dot(P.μ[i], P.Eq[i])
  end


  for mol in p.mols
    c, hs = mol[1], mol[2:5]

    for i in hs
      E += _Morse!(forces, u, c, i, P.D, P.a, P.req)
    end

    for i in 1:4
      for j in i+1:4
        E += _harmonicBondAngle!(forces, u, hs[i], c, hs[j], P.K, P.θeq)
      end
    end
    
  end

  for par in p.pars
    c1, hs1 = par[1][1], par[1][2:5]
    c2, hs2 = par[2][1], par[2][2:5] 

    E += _interTTM!(
      forces, u, P.μ, c1, c2, P.Qc, P.Qc, 
      P.Acc, P.Bcc, P.Ccc, P.A6cc; damp=tangToennies, p=P.Bcc
    )

    k = true
    for i in hs1
      E += _interTTM!(
        forces, u, P.μ, i, c2, P.Qh, P.Qc, 
        P.Ach, P.Bch, P.Cch, P.A6ch; damp=tangToennies, p=P.Bch
      )
      for j in hs2
        E += _interTTM!(
          forces, u, P.μ, i, j, P.Qh, P.Qh, 
          P.Ahh, P.Bhh, P.Chh, P.A6hh; damp=tangToennies, p=P.Bhh
        )
        if k
          E += _interTTM!(
            forces, u, P.μ, j, c1, P.Qh, P.Qc, 
            P.Ach, P.Bch, P.Cch, P.A6ch; damp=tangToennies, p=P.Bch
          )
        end

      end
      k = false
    end
  end

  if G != nothing
    tmp = [j for i in forces for j in i]
    for i in 1:length(G)
      G[i] = -tmp[i]
    end
  end

  if F != nothing
    return E
  end

end