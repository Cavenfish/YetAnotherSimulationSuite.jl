"""
TIP4P/2005f 
"""

function TIP4P(dv, v, u, p, t)
  D   = 4.48339    # eV
  a   = 2.287      # \AA
  req = 0.9419     # \AA
  K   = 3.537308e-7# eV rad^-2  1.16123e-3 * (2pi/360)^2
  θeq = 1.87448    # rad        107.4  * (2pi/360)
  ϵoo = 8.03e-3    # eV
  σoo = 3.1644     # \AA
  dm  = 0.13194    # \AA
  Qh  = 2.1113635  # 
  Qm  = - 2Qh      #
  

  # initialize things
  E = 0.0
  F = zero(u)


  for mol in p.mols
    h1, h2, o     = mol

    E += _Morse!(F, u, o, h1, D, a, req)
    E += _Morse!(F, u, o, h2, D, a, req)
    E += _harmonicBondAngle!(F, u, h1, o, h2, K, θeq)
  end

  for par in p.pars
    h1, h2, o1 = par[1]
    h3, h4, o2 = par[2]

    E += _vdw!(F, u, o1, o2, ϵoo, σoo)

    for i in [h1, h2, m1]
      for j in [h3, h4, m2]
        E += _Coulomb!(F, u, i, j, Qi, Qj)
      end
    end

  end


  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end


function TIP4P(F, G, y0, p)
  D   = 4.48339    # eV
  a   = 2.287      # \AA
  req = 0.9419     # \AA
  K   = 3.81209321 # eV rad^-2  1.16123e-3 * (360/2pi)^2
  θeq = 1.87448    # rad        107.4  * (2pi/360)
  ϵoo = 8.03e-3    # eV
  σoo = 3.1644     # \AA
  dm  = 0.13194    # \AA
  Qh  = 2.1113635  # 
  Qm  = - 2Qh      #
  
  # initialize things
  E      = 0.0
  u      = Vector[]
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

  for mol in p.mols
    h1, h2, o     = mol

    E += _Morse!(forces, u, o, h1, D, a, req)
    E += _Morse!(forces, u, o, h2, D, a, req)
    E += _harmonicBondAngle!(forces, u, h1, o, h2, K, θeq)
  end

  for par in p.pars
    h1, h2, o1 = par[1]
    h3, h4, o2 = par[2]

    E += _vdw!(forces, u, o1, o2, ϵoo, σoo)

    for i in [h1, h2, m1]
      for j in [h3, h4, m2]
        E += _Coulomb!(forces, u, i, j, Qi, Qj)
      end
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