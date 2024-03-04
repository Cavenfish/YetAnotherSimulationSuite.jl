"""
TIP4P/2005f 
"""

function TIP4P(dv, v, u, p, t)
  D   =
  a   =
  req =
  K   =
  θeq = 
  

  # initialize things
  E = 0.0
  F = zero(u)


  for mol in p.mols
    h1, h2, o     = u[mol]

    E += _Morse(F, u, o, h1, D, a, req)
    E += _Morse(F, u, o, h2, D, a, req)

    e3, F1, F2, Fo = _harmonicBondAngle()
  end



  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

end