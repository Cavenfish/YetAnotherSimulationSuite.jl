"""
Intermolecular Potential Functions
"""


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
  F[i] -= f
  F[j] += f

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

function _shortDisp(ri, rj, Aij, Bij)
  rvec = rj - ri
  r    = norm(rvec)
  E    = Aij * exp(-Bij * r)
  F    = Bij * E * rvec / r

  E,F
end

function _shortDisp!(F, u, i, j, Aij, Bij)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = Aij * exp(-Bij * r)
  f     = Bij * E * rvec / r
  F[i] -= f
  F[j] += f

  E
end

function _longDisp(ri, rj, Cij; damp=nothing, p=nothing)
  rvec = rj - ri
  r    = norm(rvec)
  
  # Get damping value if wanted
  damp != nothing ? d = damp(r, p) : d = 1 

  E = Cij / r^6 * d
  F = 6 * E * rvec / r^2

  E,F
end

function _longDisp!(F, u, i, j, Cij; damp=nothing, p=nothing)
  rvec  = u[j] - u[i]
  r     = norm(rvec)

  # Get damping value if wanted
  damp != nothing ? d = damp(r, p) : d = 1 

  E     = Cij / r^6 * d
  f     = 6 * E * rvec / r^2
  F[i] -= f
  F[j] += f

  E
end
