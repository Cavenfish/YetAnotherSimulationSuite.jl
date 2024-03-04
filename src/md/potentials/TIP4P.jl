"""
TIP4P/2005f 
"""

function _harmonicBondAngle(r1, r2, K, θeq)
  θ   = dot(r1, r2) / (norm(r1) * norm(r2)) |> acos
  E   = 0.5 * K * (θ - θeq)^2
  pre = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(r1) * norm(r2))
  F1  = pre * (r2 - (r1 * (dot(r1, r2) / dot(r1,r1))))
  F2  = pre * (r1 - (r2 * (dot(r1, r2) / dot(r2,r2))))
  Fo  = - (F1 + F2)

  E, F1, F2, Fo
end
  

function TIP4P(dv, v, u, p, t)




  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

end