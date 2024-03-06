
function _vdw(ri, rj, ϵij, σij)
  rij = norm(rj - ri)
  a   = σij / (rij)
  E   = 4ϵij * ((a)^12 - (a)^6)
  F   = 4ϵij * (12*(a)^11 - 6*(a)^5) * (σij / rij^3)

  E,F
end

function _vdw!(F, u, i, j, ϵij, σij)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  a     = σij / r
  E     = 4ϵij * ((a)^12 - (a)^6)
  f     = 4ϵij * (12*(a)^11 - 6*(a)^5) * (σij / r^3) * rvec
  F[i] += f
  F[j] -= f

  E
end

function _Coulomb(ri, rj, Qi, Qj)
  rvec = rj - ri
  r    = norm(rvec)
  E    = Qi*Qj / r
  F    = Qi*Qj * rvec / r^3

  E,F
end

function _Coulomb!(F, u, i, j, Qi, Qj)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = Qi*Qj / r
  f     = Qi*Qj * rvec / r^3
  F[i] -= f
  F[j] += f

  E
end