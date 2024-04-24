"""
Intermolecular Potential Functions
"""


function _vdw(ri, rj, ϵij, σij)
  rvec = rj - ri
  r    = norm(rvec)
  a    = σij / (r)
  E    = 4ϵij * ((a)^12 - (a)^6)
  F    = 4ϵij * (12*(a)^11 - 6*(a)^5) * (σij / r^3) * rvec

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

function _Buckingham(ri, rj, Aij, Bij, Cij)
  rvec = rj - ri
  r    = norm(rvec)
  a    = Aij * exp(-Bij * r)
  b    = Cij / r^6
  E    = a - b
  F    = (Bij * a / r * rvec) - (6b / r^2 * rvec)

  E,F
end

function _Buckingham!(F, u, i, j, Aij, Bij, Cij)
  rvec = u[j] - u[i]
  r    = norm(rvec)
  a    = Aij * exp(-Bij * r)
  b    = Cij / r^6
  E    = a - b
  f    = (Bij * a / r * rvec) - (6b / r^2 * rvec)

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

  E = -Cij / r^6 * d
  F = 6 * E * rvec / r^2

  E,F
end

function _longDisp!(F, u, i, j, Cij; damp=nothing, p=nothing)
  rvec  = u[j] - u[i]
  r     = norm(rvec)

  # Get damping value if wanted
  damp != nothing ? d = damp(r, p) : d = 1
  
  E     = -Cij / r^6 * d
  f     = 6 * E * rvec / r^2
  F[i] -= f
  F[j] += f

  E
end

function _Vpol4Fcc(ri, rj, Qi, Qj, A)
  rvec = rj - ri
  r    = norm(rvec)
  a    = 0.4

  c    = a * (r/A)^4
  s0   = 1 - exp(-c) + (a^0.25 * r / A) * gamma(0.75, c)
  s1   = 1 - exp(-c)

  Eqq  = s0 * Qi * Qj / r
  Fqq  = s1 * Qi * Qj / r^3 * rvec

  Eqq, Fqq
end

function _Vpol4Fcc!(F, u, i, j, Qi, Qj, A)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  a     = 0.4

  c     = a * (r/A)^4
  s0    = 1 - exp(-c) + (a^0.25 * r / A) * gamma(0.75, c)
  s1    = 1 - exp(-c)

  Eqq   = s0 * Qi * Qj / r
  Fqq   = s1 * Qi * Qj / r^3 * rvec

  F[i] -= Fqq
  F[j] += Fqq

  Eqq
end

function _Vpol4Fcd(ri, rj, Qi, Qj, μi, μj, A)
  rvec = rj - ri
  r    = norm(rvec)
  a    = 0.4

  c     = exp(-a*(r/A)^4)
  s1    =  1 - c
  s2    = s1 - (4a/3) * (r/A)^4 * c
  s3    = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi   = dot(rvec, μi)
  rμj   = dot(rvec, μj)
  Equ   = (Qi * rμj - Qj * rμi) * (s1 / r^3)
  Fqu   = (Qi * rμj - Qj * rμi) * (s2 * 3 / r^5 * rvec) .+ (Qj * μi .- Qi * μj) * (s1 / r^3)

  Equ, Fqu
end

function _Vpol4Fcd!(F, u, i, j, Qi, Qj, μi, μj, A)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  a     = 0.4

  c     = exp(-a*(r/A)^4)
  s1    =  1 - c
  s2    = s1 - (4a/3) * (r/A)^4 * c
  s3    = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi   = dot(rvec, μi)
  rμj   = dot(rvec, μj)
  Equ   = (Qi * rμj - Qj * rμi) * (s1 / r^3)
  Fqu   = (Qi * rμj - Qj * rμi) * (s2 * 3 / r^5 * rvec) .+ (Qj * μi .- Qi * μj) * (s1 / r^3)

  F[i] -= Fqu
  # F[j] += Fqu

  Equ
end

function _Vpol4Fdd(ri, rj, Qi, Qj, μi, μj, A)
  rvec = rj - ri
  r    = norm(rvec)
  a    = 0.055*0.626

  c    = exp(-a*(r/A)^4)
  s1   =  1 - c
  s2   = s1 - (4a/3) * (r/A)^4 * c
  s3   = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi  = dot(rvec, μi)
  rμj  = dot(rvec, μj)
  μij  = dot(μi, μj)

  Euu  = (s1 / r^3 * μij) - (s2 * 3 / r^5 * rμi * rμj)
  Fuu  = (-s3 * 15 / r^7 * rμi * rμj * rvec) .+ (s2 * 3 / r^5) * (μij * rvec .+ rμi * μj .+ rμj * μi)

  Euu, Fuu
end

function _Vpol4Fdd!(F, u, i, j, Qi, Qj, μi, μj, A)
  rvec = u[j] - u[i]
  r    = norm(rvec)
  a    = 0.055*0.626

  c    = exp(-a*(r/A)^4)
  s1   =  1 - c
  s2   = s1 - (4a/3) * (r/A)^4 * c
  s3   = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi  = dot(rvec, μi)
  rμj  = dot(rvec, μj)
  μij  = dot(μi, μj)

  Euu  = (s1 / r^3 * μij) - (s2 * 3 / r^5 * rμi * rμj)
  Fuu  = (-s3 * 15 / r^7 * rμi * rμj * rvec) .+ (s2 * 3 / r^5) * (μij * rvec .+ rμi * μj .+ rμj * μi)

  F[i] -= Fuu
  # F[j] += Fuu

  Euu
end