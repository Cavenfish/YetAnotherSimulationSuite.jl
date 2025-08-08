"""
Intramolecular Potential Functions
"""


function _Morse(
  r::Float64, rvec::V, D::Float64, a::Float64, req::Float64
) where V <: AbstractVector

  c = exp(-a*(r-req))
  E = D * (1 - c)^2
  F = @. -2D * a * c * (1 - c) * rvec / r

  E, F
end

function _Morse!(
  F::Vector{Vf}, u::Vector{Vu}, 
  i::Int64, j::Int64, D::Float64, a::Float64, req::Float64
) where {Vf <: AbstractVector, Vu <: AbstractVector}
  
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  c     = exp(-a*(r-req))
  E     = D * (1 - c)^2
  f     = @. -2D * a * c * (1 - c) * rvec / r

  F[i] .-= f
  F[j] .+= f

  E
end

function _Morse!(
  F::Vector{Vf}, u::Vector{Vu}, Fbuf::BUF, rbuf::BUF,
  i::Int64, j::Int64, D::Float64, a::Float64, req::Float64
) where {Vf <: AbstractVector, Vu <: AbstractVector, BUF<:AbstractVector}
  
  @. rbuf  = u[j] - u[i]
  r     = norm(rbuf)
  c     = exp(-a*(r-req))
  E     = D * (1 - c)^2
  @. Fbuf  = -2D * a * c * (1 - c) * rbuf / r

  F[i] .-= Fbuf
  F[j] .+= Fbuf

  E
end

function _harmonicBond(
  r::Float64, rvec::V, K::Float64, req::Float64
) where V <: AbstractVector

  E = 0.5 * K * (r - req)^2
  F = @. - K * (r - req) * rvec / r

  E, F
end

function _harmonicBond!(
  F::Vector{Vf}, u::Vector{Vu}, 
  i::Int64, j::Int64, K::Float64, req::Float64
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = 0.5 * K * (r - req)^2
  f     = @. - K * (r - req) * rvec / r

  F[i] .-= f
  F[j] .+= f

  E
end

function _harmonicBond!(
  F::Vector{Vf}, u::Vector{Vu}, Fbuf::BUF, rbuf::BUF,
  i::Int64, j::Int64, K::Float64, req::Float64
) where {Vf <: AbstractVector, Vu <: AbstractVector, BUF<:AbstractVector}

  @. rbuf  = u[j] - u[i]
  r     = norm(rbuf)
  E     = 0.5 * K * (r - req)^2
  @. Fbuf = -K * (r - req) * rbuf / r

  F[i] .-= Fbuf
  F[j] .+= Fbuf

  E
end

function _harmonicBondAngle(
  r1::V, r2::V, K::Float64, θeq::Float64
) where V <: AbstractVector

  θ   = dot(r1, r2) / (norm(r1) * norm(r2)) |> (x -> round(x, digits=10)) |> acos
  E   = 0.5 * K * (θ - θeq)^2
  pre = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(r1) * norm(r2))
  F1  = pre * (r2 - (r1 * (dot(r1, r2) / dot(r1,r1))))
  F2  = pre * (r1 - (r2 * (dot(r1, r2) / dot(r2,r2))))
  Fo  = @. - (F1 + F2)

  E, F1, F2, Fo
end

function _harmonicBondAngle!(
  F::Vector{Vf}, u::Vector{Vu}, 
  i::Int64, o::Int64, j::Int64, K::Float64, θeq::Float64
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  ri    = u[i] - u[o]
  rj    = u[j] - u[o]
  θ     = dot(ri, rj) / (norm(ri) * norm(rj)) |> (x -> round(x, digits=10)) |> acos
  E     = 0.5 * K * (θ - θeq)^2
  pre   = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(ri) * norm(rj))
  Fi    = pre * (rj - (ri * (dot(ri, rj) / dot(ri,ri))))
  Fj    = pre * (ri - (rj * (dot(ri, rj) / dot(rj,rj))))
  F[i] .+= Fi
  F[j] .+= Fj
  @. F[o] -= (Fi + Fj)

  E
end