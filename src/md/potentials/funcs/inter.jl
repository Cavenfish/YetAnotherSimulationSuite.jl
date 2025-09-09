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
  F::Vector{Vf}, u::Vector{Vu},
  i::Int64, j::Int64, ϵ::Float64, σ::Float64; S=1.0
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  rvec  = u[j] - u[i]
  r     = norm(rvec)
  a     = σ / r
  E     = 4ϵ * ((a)^12 - (a)^6)
  f     = 4ϵ * (12*(a)^11 - 6*(a)^5) * (σ / r^3) * rvec
  @. F[i] -= f * S
  @. F[j] += f * S

  E
end

function _vdw!(
  F::Vector{Vf}, u::Vector{Vu}, Fbuf::BUF, rbuf::BUF,
  i::Int64, j::Int64, ϵ::Float64, σ::Float64; S=1.0
) where {Vf <: AbstractVector, Vu <: AbstractVector, BUF<:AbstractVector}

  @. rbuf = u[j] - u[i]
  r     = norm(rbuf)
  a     = σ / r
  E     = 4ϵ * ((a)^12 - (a)^6)
  @. Fbuf  = 4ϵ * (12*(a)^11 - 6*(a)^5) * (σ / r^3) * rbuf
  @. F[i] -= Fbuf * S
  @. F[j] += Fbuf * S

  E
end

function _vdw!(
  F::Vector{Vf}, u::Vector{Vu}, Fbuf::BUF, rbuf::BUF, lat::AbstractMatrix,
  i::Int64, j::Int64, ϵ::Float64, σ::Float64, rc::Float64; S=1.0
) where {Vf <: AbstractVector, Vu <: AbstractVector, BUF<:AbstractVector}

  r        = pbcVec!(rbuf, u[i], u[j], rc, lat)
  a        = σ / r
  E        = 4ϵ * ((a)^12 - (a)^6)
  @. Fbuf  = 4ϵ * (12*(a)^11 - 6*(a)^5) * (σ / r^3) * rbuf
  @. F[i] -= Fbuf * S
  @. F[j] += Fbuf * S

  E
end

function _Buckingham(
  ri::Vi, rj::Vj, Aij::Float64, Bij::Float64, Cij::Float64
) where {Vi <: AbstractVector, Vj <: AbstractVector}

  rvec = rj - ri
  r    = norm(rvec)
  a    = Aij * exp(-Bij * r)
  b    = Cij / r^6
  E    = a - b
  F    = @. (Bij * a / r * rvec) - (6b / r^2 * rvec)

  E,F
end

function _Buckingham!(
  F::Vector{Vf}, u::Vector{Vu},
  i::Int64, j::Int64, Aij::Float64, Bij::Float64, Cij::Float64
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  rvec = u[j] - u[i]
  r    = norm(rvec)
  a    = Aij * exp(-Bij * r)
  b    = Cij / r^6
  E    = a - b
  f    = @. (Bij * a / r * rvec) - (6b / r^2 * rvec)

  F[i] .-= f
  F[j] .+= f

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
  F::BUF, rbuf::BUF, ri::Vi, rj::Vj, Qi::Float64, Qj::Float64
) where {Vi <: AbstractVector, Vj <: AbstractVector, BUF<:AbstractVector}

  rbuf = rj - ri
  r    = norm(rbuf)
  E    = Qi*Qj / r
  @. F = Qi*Qj / r^3 * rbuf

  E
end

function _Coulomb!(
  F::Vector{Vf}, u::Vector{Vu}, 
  i::Int64, j::Int64, Qi::Float64, Qj::Float64; S=1.0
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = Qi*Qj / r
  f     = @. Qi*Qj / r^3 * rvec
  @. F[i] .-= f * S
  @. F[j] .+= f * S

  E
end

function _Coulomb!(
  F::Vector{Vf}, u::Vector{Vu}, Fbuf::BUF, rbuf::BUF,
  i::Int64, j::Int64, Qi::Float64, Qj::Float64; S=1.0
) where {Vf <: AbstractVector, Vu <: AbstractVector, BUF<:AbstractVector}

  @. rbuf   = u[j] - u[i]
  r         = norm(rbuf)
  E         = Qi*Qj / r
  @. Fbuf   = Qi*Qj / r^3 * rbuf
  @. F[i] .-= Fbuf * S
  @. F[j] .+= Fbuf * S

  E
end

function _Coulomb!(
  F::Vector{Vf}, u::Vector{Vu}, Fbuf::BUF, rbuf::BUF, lat::AbstractMatrix,
  i::Int64, j::Int64, Qi::Float64, Qj::Float64, rc::Float64; S=1.0
) where {Vf <: AbstractVector, Vu <: AbstractVector, BUF<:AbstractVector}

  r         = pbcVec!(rbuf, u[i], u[j], rc, lat)
  E         = Qi*Qj / r
  @. Fbuf   = Qi*Qj / r^3 * rbuf
  @. F[i] .-= Fbuf * S
  @. F[j] .+= Fbuf * S

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
  F::Vector{Vf}, u::Vector{Vu},
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
  F::Vector{Vf}, u::Vector{Vu}, 
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
