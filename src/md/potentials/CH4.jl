"""
CH4 MBnrg Potential
from Riera et al. 2020
https://doi.org/10.1021/acs.jpcb.0c08728
"""


function CH4(dv, v, u, p, t)
  D   = 4.0017 # eV
  a   = 2.0415 # \AA^-1
  req = 1.0887 # \AA
  K   = 3.3349 # eV rad^-1
  θeq = 1.8741 # rad 107.379 * (2pi/360)

  # initialize things
  E = 0.0
  F = zero(u)


  for mol in p.mols
    h1, h2, h3, h4, c = mol

    for i in [h1, h2, h3, h4]
      E += _Morse!(F, u, c, i, D, a, req)
    end

    for i in [h1, h2, h3, h4]
      for j in [h1, h2, h3, h4]
        i == j && continue 
        E += _harmonicBondAngle!(F, u, i, c, j, K, θeq)
      end
    end
    
  end

  for par in p.pars
  end
  

  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end

function CH4(F, G, y0, p)
  D   = 4.0017 # eV
  a   = 2.0415 # \AA^-1
  req = 1.0887 # \AA
  K   = 3.3349 # eV rad^-1
  θeq = 1.8741 # rad 107.379 * (2pi/360)

  # initialize things
  E      = 0.0
  u      = Vector[]
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

  for mol in p.mols
    c  = mol[1]
    hs = mol[2:5]

    for i in hs
      E += _Morse!(forces, u, c, i, D, a, req) - D
    end

    for i in 1:4
      for j in i+1:4
        E += _harmonicBondAngle!(forces, u, hs[i], c, hs[j], K, θeq) - K
      end
    end
    
  end

  for par in p.pars
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