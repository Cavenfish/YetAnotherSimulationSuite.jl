"""
CH4 TTM-nrg Potential
from Riera et al. 2020
https://doi.org/10.1021/acs.jpcb.0c08728
"""


function CH4(dv, v, u, p, t)
  D   = 4.0017  # eV
  a   = 2.0415  # \AA^-1
  req = 1.0887  # \AA
  K   = 3.3349  # eV rad^-1
  θeq = 1.8741  # rad 107.379 * (2pi/360)
  Acc = 1852.25 #
  Ach = 141.317 #
  Ahh = 112.503 #
  Bcc = 3.37925 # \AA^-1
  Bch = 3.25885 # \AA^-1
  Bhh = 4.05972 # \AA^-1
  Ccc = 13.15   #
  Cch = 4.51455 #
  Chh = 1.59497 #
  Qh  = 0.51163 #
  Qc  = - 4Qh   #
  αh  = 0.38978 # \AA^3
  αc  = 1.39327 # \AA^3

  # initialize things
  E = 0.0
  F = zero(u)


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

    E += _Coulomb!(F, u, c1, c2, Qc, Qc)
    E += _shortDisp!(F, u, c1, c2, Acc, Bcc)
    E += _longDisp!(F, u, c1, c2, Ccc; damp=tangToennies, p=Bcc)

    k = true
    for i in hs1
      E += _Coulomb!(F, u, i, c2, Qh, Qc)
      E += _shortDisp!(F, u, i, c2, Ach, Bch)
      E += _longDisp!(F, u, i, c2, Cch; damp=tangToennies, p=Bch)
      for j in hs2
        E += _Coulomb!(F, u, i, j, Qh, Qh)
        E += _shortDisp!(F, u, i, j, Ahh, Bhh)
        E += _longDisp!(F, u, i, j, Chh; damp=tangToennies, p=Bhh)
        
        if k
          E += _Coulomb!(F, u, j, c1, Qh, Qc)
          E += _shortDisp!(F, u, j, c1, Ach, Bch)
          E += _longDisp!(F, u, j, c1, Cch; damp=tangToennies, p=Bch)
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
  D   = 4.0017  # eV
  a   = 2.0415  # \AA^-1
  req = 1.0887  # \AA
  K   = 3.3349  # eV rad^-1
  θeq = 1.8741  # rad 107.379 * (2pi/360)
  Acc = 1852.25 #
  Ach = 141.317 #
  Ahh = 112.503 #
  Bcc = 3.37925 # \AA^-1
  Bch = 3.25885 # \AA^-1
  Bhh = 4.05972 # \AA^-1
  Ccc = 13.15   #
  Cch = 4.51455 #
  Chh = 1.59497 #
  Qh  = 0.51163 #
  Qc  = - 4Qh   #
  αh  = 0.38978 # \AA^3
  αc  = 1.39327 # \AA^3

  # initialize things
  E      = 0.0
  u      = Vector[]
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

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

    E += _Coulomb!(forces, u, c1, c2, Qc, Qc)
    E += _shortDisp!(forces, u, c1, c2, Acc, Bcc)
    E += _longDisp!(forces, u, c1, c2, Ccc; damp=tangToennies, p=Bcc)

    E += _Vpol4Fcc!(forces, u, c1, c2, Qc, Qc, αc^(2/6))

    k = true
    for i in hs1
      E += _Coulomb!(forces, u, i, c2, Qh, Qc)
      E += _shortDisp!(forces, u, i, c2, Ach, Bch)
      E += _longDisp!(forces, u, i, c2, Cch; damp=tangToennies, p=Bch)
      E += _Vpol4Fcc!(forces, u, i, c2, Qh, Qc, (αc*αh)^(1/6))
      for j in hs2
        E += _Coulomb!(forces, u, i, j, Qh, Qh)
        E += _shortDisp!(forces, u, i, j, Ahh, Bhh)
        E += _longDisp!(forces, u, i, j, Chh; damp=tangToennies, p=Bhh)
        E += _Vpol4Fcc!(forces, u, i, j, Qh, Qh, αh^(2/6))
        
        if k
          E += _Coulomb!(forces, u, j, c1, Qh, Qc)
          E += _shortDisp!(forces, u, j, c1, Ach, Bch)
          E += _longDisp!(forces, u, j, c1, Cch; damp=tangToennies, p=Bch)
          E += _Vpol4Fcc!(forces, u, j, c1, Qh, Qc, (αc*αh)^(1/6))
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

function _getDipole(mol, u, m, q)
  o = CoM(u[mol], m[mol])
  μ = zeros(3)
  
  for i in 1:length(mol)
    j  = mol[i]
    r  = u[j] - o
    μ += q[i] * r
  end

  μ
end

