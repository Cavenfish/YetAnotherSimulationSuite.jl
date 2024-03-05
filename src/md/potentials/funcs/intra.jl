

function _Morse(r, rvec, D, a, req)
  preF = 2 * D * a / req
  expf = exp(a * (1 - r / req))
  E    = D * expf * (expf - 2)
  F    = preF * expf * (expf - 1) * rvec / r

  E, F
end

function _Morse!(F, u, i, j, D, a, req)
  r, rvec = diffDotSqrt(u[j], u[i])
  preF    = 2 * D * a / req
  expf    = exp(a * (1 - r / req))
  E       = D * expf * (expf - 2)
  f       = preF * expf * (expf - 1) * rvec / r

  F[i] -= f
  F[j] += f

  E
end

function _harmonicBondAngle(r1, r2, K, θeq)
  θ   = dot(r1, r2) / (norm(r1) * norm(r2)) |> acos
  E   = 0.5 * K * (θ - θeq)^2
  pre = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(r1) * norm(r2))
  F1  = pre * (r2 - (r1 * (dot(r1, r2) / dot(r1,r1))))
  F2  = pre * (r1 - (r2 * (dot(r1, r2) / dot(r2,r2))))
  Fo  = - (F1 + F2)

  E, F1, F2, Fo
end

function _harmonicBondAngle!(F, u, i, o, j, K, θeq)
  ri    = u[i]
  rj    = u[j]
  θ     = dot(ri, rj) / (norm(ri) * norm(rj)) |> acos
  E     = 0.5 * K * (θ - θeq)^2
  pre   = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(ri) * norm(rj))
  F[i] += pre * (rj - (ri * (dot(ri, rj) / dot(ri,ri))))
  F[j] += pre * (ri - (rj * (dot(ri, rj) / dot(rj,rj))))
  F[o] -= (Fi + Fj)

  E
end