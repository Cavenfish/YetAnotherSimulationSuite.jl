"""
CH4 TTM-nrg Potential
from Riera et al. 2020
https://doi.org/10.1021/acs.jpcb.0c08728
"""

struct _CH4_PotVars <: PotVars
  μ::Vector
end

function CH4(bdys::Vector{Atom})
  μ = [zeros(3) for i in bdys]

  _CH4_PotVars(μ)
end

function CH4(dv, v, u, p, t)
  D    = 4.0017  # eV
  a    = 2.0415  # \AA^-1
  req  = 1.0887  # \AA
  K    = 3.33487 # eV rad^-1   76.9061 * 0.043363
  θeq  = 1.87411 # rad 107.379 * (2pi/360)
  Acc  = 1852.25 #
  Ach  = 141.317 #
  Ahh  = 112.503 #
  Bcc  = 3.37925 # \AA^-1
  Bch  = 3.25885 # \AA^-1
  Bhh  = 4.05972 # \AA^-1
  Ccc  = 13.15   #
  Cch  = 4.51455 #
  Chh  = 1.59497 #
  Qh   = 0.51163 #
  Qc   = - 4Qh   #
  αh   = 0.38978 # \AA^3
  αc   = 1.39327 # \AA^3
  A6cc = αc^(2/6)
  A6ch = (αc*αh)^(1/6)
  A6hh = αh^(2/6)

  # initialize things
  E = 0.0
  F = zero(u)
  α = repeat([αc, αh, αh, αh, αh], length(p.mols))
  Q = repeat([Qc, Qh, Qh, Qh, Qh], length(p.mols))

  _getDipoles4TTM_MatrixInversion!(p.μ, u, Q, α)

  E += _getDipolePolarizationEnergy(p.μ, α)

  for mol in p.mols
    c, hs = mol[1], mol[2:5]

    for i in hs
      E += _Morse!(F, u, c, i, D, a, req)
    end

    for i in 1:4
      for j in i+1:4
        E += _harmonicBondAngle!(F, u, hs[i], c, hs[j], K, θeq)
      end
    end
    
  end

  for par in p.pars
    c1, hs1 = par[1][1], par[1][2:5]
    c2, hs2 = par[2][1], par[2][2:5] 

    E += _interTTM!(
      F, u, p.μ, c1, c2, Qc, Qc, 
      Acc, Bcc, Ccc, A6cc; damp=tangToennies, p=Bcc
    )

    k = true
    for i in hs1
      E += _interTTM!(
        F, u, p.μ, i, c2, Qh, Qc, 
        Ach, Bch, Cch, A6ch; damp=tangToennies, p=Bch
      )
      for j in hs2
        E += _interTTM!(
          F, u, p.μ, i, j, Qh, Qh, 
          Ahh, Bhh, Chh, A6hh; damp=tangToennies, p=Bhh
        )
        if k
          E += _interTTM!(
            F, u, p.μ, j, c1, Qh, Qc, 
            Ach, Bch, Cch, A6ch; damp=tangToennies, p=Bch
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
  D    = 4.0017  # eV
  a    = 2.0415  # \AA^-1
  req  = 1.0887  # \AA
  K    = 3.33487 # eV rad^-1   76.9061 * 0.043363
  θeq  = 1.87411 # rad 107.379 * (2pi/360)
  Acc  = 1852.25 #
  Ach  = 141.317 #
  Ahh  = 112.503 #
  Bcc  = 3.37925 # \AA^-1
  Bch  = 3.25885 # \AA^-1
  Bhh  = 4.05972 # \AA^-1
  Ccc  = 13.15   #
  Cch  = 4.51455 #
  Chh  = 1.59497 #
  Qh   = 0.51163 #
  Qc   = - 4Qh   #
  αh   = 0.38978 # \AA^3
  αc   = 1.39327 # \AA^3
  A6cc = αc^(2/6)
  A6ch = (αc*αh)^(1/6)
  A6hh = αh^(2/6)

  # initialize things
  α      = repeat([αc, αh, αh, αh, αh], length(p.mols))
  Q      = repeat([Qc, Qh, Qh, Qh, Qh], length(p.mols))
  E      = 0.0
  u      = Vector[]
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

  # _getDipoles4TTM_Iterative!(p.μ, u, Q, α, p.mols)
  _getDipoles4TTM_MatrixInversion!(p.μ, u, Q, α)

  E += _getDipolePolarizationEnergy(p.μ, α)

  for mol in p.mols
    c, hs = mol[1], mol[2:5]

    for i in hs
      E += _Morse!(forces, u, c, i, D, a, req)
    end

    for i in 1:4
      for j in i+1:4
        E += _harmonicBondAngle!(forces, u, hs[i], c, hs[j], K, θeq)
      end
    end
    
  end

  for par in p.pars
    c1, hs1 = par[1][1], par[1][2:5]
    c2, hs2 = par[2][1], par[2][2:5] 

    E += _interTTM!(
      forces, u, p.μ, c1, c2, Qc, Qc, 
      Acc, Bcc, Ccc, A6cc; damp=tangToennies, p=Bcc
    )

    k = true
    for i in hs1
      E += _interTTM!(
        forces, u, p.μ, i, c2, Qh, Qc, 
        Ach, Bch, Cch, A6ch; damp=tangToennies, p=Bch
      )
      for j in hs2
        E += _interTTM!(
          forces, u, p.μ, i, j, Qh, Qh, 
          Ahh, Bhh, Chh, A6hh; damp=tangToennies, p=Bhh
        )
        if k
          E += _interTTM!(
            forces, u, p.μ, j, c1, Qh, Qc, 
            Ach, Bch, Cch, A6ch; damp=tangToennies, p=Bch
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