"""
Intermolecular Potential Functions
"""


function _vdw(
  ri::Vi, rj::Vj, ϵ::Float64, σ::Float64
) where {Vi <: AbstractVector, Vj <: AbstractVector}

  rvec = rj - ri
  r    = norm(rvec)
  a    = σ / r
  E    = 4ϵ * ((a)^12 - (a)^6)
  F    = @. 4ϵ * (12*(a)^11 - 6*(a)^5) * (σ / r^3) * rvec

  E,F
end

function _vdw!(
  F::Vf, u::Vu, i::Int64, j::Int64, ϵ::Float64, σ::Float64; 
  S=1.0, buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  @. buf.ri = u[j] - u[i]
  r     = norm(buf.ri)
  a     = σ / r
  E     = 4ϵ * ((a)^12 - (a)^6)
  @. buf.fi  = 4ϵ * (12*(a)^11 - 6*(a)^5) * (σ / r^3) * buf.ri
  @. F[i] -= buf.fi * S
  @. F[j] += buf.fi * S

  E
end

function _vdw!(
  F::Vf, u::Vu, lat::AbstractMatrix, i::Int64, j::Int64,
  ϵ::Float64, σ::Float64, rc::Float64; 
  S=1.0, buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  r         = pbcVec!(buf.ri, u[i], u[j], rc, lat)
  a         = σ / r
  E         = 4ϵ * ((a)^12 - (a)^6)
  @. buf.fi = 4ϵ * (12*(a)^11 - 6*(a)^5) * (σ / r^3) * buf.ri
  @. F[i]  -= buf.fi * S
  @. F[j]  += buf.fi * S

  E
end

function _Buckingham(
  ri::Vi, rj::Vj, A::Float64, B::Float64, C::Float64
) where {Vi <: AbstractVector, Vj <: AbstractVector}

  rvec = rj - ri
  r    = norm(rvec)
  a    = A * exp(-B * r)
  b    = C / r^6
  E    = a - b
  F    = @. (B * a / r * rvec) - (6b / r^2 * rvec)

  E,F
end

function _Buckingham!(
  F::Vf, u::Vu, i::Int64, j::Int64,
  A::Float64, B::Float64, C::Float64;
  S=1.0, buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  @. buf.ri = u[j] - u[i]
  r     = norm(buf.ri)
  a    = A * exp(-B * r)
  b    = C / r^6
  E    = a - b
  f    = @. (B * a / r * buf.ri) - (6b / r^2 * buf.ri)

  F[i] .-= f * S
  F[j] .+= f * S

  E
end

function _Buckingham!(
  F::Vf, u::Vu, lat::AbstractMatrix, i::Int64, j::Int64,
  A::Float64, B::Float64, C::Float64, rc::Float64;
  S=1.0, buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  r = pbcVec!(buf.ri, u[i], u[j], rc, lat)
  a = A * exp(-B * r)
  b = C / r^6
  E = a - b
  
  @. buf.fi = (B * a / r * buf.ri) - (6b / r^2 * buf.ri)

  F[i] .-= buf.fi * S
  F[j] .+= buf.fi * S

  E
end

function _Coulomb(
  ri::Vi, rj::Vj, Qi::Float64, Qj::Float64
) where {Vi <: AbstractVector, Vj <: AbstractVector}

  rvec = rj - ri
  r    = norm(rvec)
  E    = Qi*Qj / r
  F    = @. Qi*Qj / r^3 * rvec

  E,F
end

function _Coulomb!(
  F::AbstractVector, ri::Vi, rj::Vj, lat::AbstractMatrix,
  Qi::Float64, Qj::Float64, rc::Float64;
  buf=_func_buffs
) where {Vi <: AbstractVector, Vj <: AbstractVector}

  r    = pbcVec!(buf.ri, ri, rj, rc, lat)
  E    = Qi*Qj / r
  @. F = Qi*Qj / r^3 * buf.ri

  E
end

function _Coulomb!(
  F::Vf, u::Vu, i::Int64, j::Int64, Qi::Float64, Qj::Float64; 
  S=1.0, buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  @. buf.ri = u[j] - u[i]
  r     = norm(buf.ri)
  E     = Qi*Qj / r
  @. buf.fi = Qi*Qj / r^3 * buf.ri
  @. F[i] .-= buf.fi * S
  @. F[j] .+= buf.fi * S

  E
end

function _Coulomb!(
  F::Vf, u::Vu,  lat::AbstractMatrix, i::Int64, j::Int64,
  Qi::Float64, Qj::Float64, rc::Float64; 
  S=1.0, buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  r         = pbcVec!(buf.ri, u[i], u[j], rc, lat)
  E         = Qi*Qj / r
  @. buf.fi = Qi*Qj / r^3 * buf.ri
  @. F[i] .-= buf.fi * S
  @. F[j] .+= buf.fi * S

  E
end

function _shortDisp(
  ri::Vi, rj::Vj, Aij::Float64, Bij::Float64
) where {Vi <: AbstractVector, Vj <: AbstractVector}

  rvec = rj - ri
  r    = norm(rvec)
  E    = Aij * exp(-Bij * r)
  F    = @. Bij * E * rvec / r

  E,F
end

function _shortDisp!(
  F::Vf, u::Vu,
  i::Int64, j::Int64, Aij::Float64, Bij::Float64
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = Aij * exp(-Bij * r)
  f     = @. Bij * E * rvec / r
  F[i] .-= f
  F[j] .+= f

  E
end

function _longDisp(
  ri::Vi, rj::Vj, Cij::Float64; damp=nothing, p=nothing
) where {Vi <: AbstractVector, Vj <: AbstractVector}

  rvec = rj - ri
  r    = norm(rvec)
  
  # Get damping value if wanted
  damp != nothing ? d = damp(r, p) : d = (1,0) 

  a = -Cij / r^6
  b = -6 * a * rvec / r^2
  E = a * d[1]
  F = @. -(d[2] * a * rvec/r .+ d[1] * b)

  E,F
end

function _longDisp!(
  F::Vf, u::Vu, 
  i::Int64, j::Int64, Cij::Float64; damp=nothing, p=nothing
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  rvec  = u[j] - u[i]
  r     = norm(rvec)

  # Get damping value if wanted
  damp != nothing ? d = damp(r, p) : d = (1,0)

  a     = -Cij / r^6
  b     = -6 * a * rvec / r^2
  E     = a * d[1]
  f     = @. -(d[2] * a * rvec/r .+ d[1] * b)
  F[i] .-= f
  F[j] .+= f

  E
end
