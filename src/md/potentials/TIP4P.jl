"""
TIP4P/2005f 
"""

function TIP4P(dv, v, u, p, t)
  D   =
  a   =
  req =
  K   =
  θeq =
  ϵoo =
  σoo =
  

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

end