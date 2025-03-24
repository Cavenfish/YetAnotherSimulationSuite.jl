function diffDotSqrt(v2, v1)
  rvec = v2 - v1
  r    = sqrt(dot(rvec, rvec))
  (r, rvec)
end

function getAngle(r1, r2)
  x  = (dot(r1, r2)) / (norm(r1) * norm(r2))
  x  = round(x, digits=15)
  
  acos(x)
end

