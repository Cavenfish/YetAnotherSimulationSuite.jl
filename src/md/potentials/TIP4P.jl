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

    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(F, u, i, j, Qh, Qh)
      end
    end

    E += _getMforces!(F, u, par[1], par[2], dm, Qh, Qm)

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

    for i in [h1, h2]
      for j in [h3, h4]
        E += _Coulomb!(forces, u, i, j, Qh, Qh)
      end
    end

    E += _getMforces!(forces, u, par[1], par[2], dm, Qh, Qm)

  end

  E += D * length(p.mols) * 2
  E += K * length(p.mols)

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

function _getMforces!(F, u, w1, w2, dm, Qh, Qm)
  h1, h2, o1 = w1
  h3, h4, o2 = w2

  # Get Norms 
  d11 = norm(u[h1] - u[o1])
  d21 = norm(u[h2] - u[o1])
  d32 = norm(u[h3] - u[o2])
  d42 = norm(u[h4] - u[o2])

  # Get weights
  wh1 = dm * 0.5 * (1 + d21/d11)
  wh2 = dm * 0.5 * (1 + d11/d21)
  wh3 = dm * 0.5 * (1 + d42/d32)
  wh4 = dm * 0.5 * (1 + d32/d42)

  # Get M positions
  m1 = u[o1] + wh1*(u[h1] - u[o1]) + wh2*(u[h2] - u[o1])
  m2 = u[o2] + wh3*(u[h3] - u[o2]) + wh4*(u[h4] - u[o2])

  # H1 -- M2
  E,f    = _Coulomb(u[h1], m2, Qh, Qm)
  F[h1] -= f
  F[h3] += f * wh3
  F[h4] += f * wh4
  F[o2] += f * (1 - wh3 - wh4)

  # H2 -- M2
  e,f    = _Coulomb(u[h2], m2, Qh, Qm)
  E     += e
  F[h2] -= f
  F[h3] += f * wh3
  F[h4] += f * wh4
  F[o2] += f * (1 - wh3 - wh4)

  # H3 -- M1
  e,f    = _Coulomb(u[h3], m1, Qh, Qm)
  E     += e
  F[h3] -= f
  F[h1] += f * wh1
  F[h2] += f * wh2
  F[o1] += f * (1 - wh1 - wh2)

  # H4 -- M1
  e,f    = _Coulomb(u[h4], m1, Qh, Qm)
  E     += e
  F[h4] -= f
  F[h1] += f * wh1
  F[h2] += f * wh2
  F[o1] += f * (1 - wh1 - wh2)

  # M1 -- M2
  e,f    = _Coulomb(m1, m2, Qm, Qm)
  E     += e
  F[h1] -= f * wh1
  F[h2] -= f * wh2
  F[o1] -= f * (1 - wh1 - wh2)
  F[h3] += f * wh3
  F[h4] += f * wh4
  F[o2] += f * (1 - wh3 - wh4)

  E
end