"""
Intramolecular Potential Functions
"""

struct _Buffers{AV3D<:AbstractVector}
  ri::AV3D
  rj::AV3D
  fi::AV3D
  fj::AV3D
end

const _func_buffs = _Buffers(
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3)),
  MVector{3}(zeros(3))
)

function _Morse(
  r::Float64, rvec::V, D::Float64, a::Float64, req::Float64
) where V <: AbstractVector

  c = exp(-a*(r-req))
  E = D * (1 - c)^2
  F = @. -2D * a * c * (1 - c) * rvec / r

  E, F
end

function _Morse!(
  F::Vf, u::Vu, i::Int64, j::Int64,
  D::Float64, a::Float64, req::Float64;
  buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}
  
  @. buf.ri  = u[j] - u[i]
  r     = norm(buf.ri)
  c     = exp(-a*(r-req))
  E     = D * (1 - c)^2
  @. buf.fi  = -2D * a * c * (1 - c) * buf.ri / r

  F[i] .-= buf.fi
  F[j] .+= buf.fi

  E
end

function _Morse!(
  F::Vf, u::Vu, lat::AbstractMatrix, i::Int64, j::Int64,
  D::Float64, a::Float64, req::Float64, rc::Float64;
  buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}
  
  r = pbcVec!(buf.ri, u[i], u[j], rc, lat)
  c = exp(-a*(r-req))
  E = D * (1 - c)^2

  @. buf.fi = -2D * a * c * (1 - c) * buf.ri / r

  F[i] .-= buf.fi
  F[j] .+= buf.fi

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
  F::Vf, u::Vu, i::Int64, j::Int64, K::Float64, req::Float64;
  buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  @. buf.ri  = u[j] - u[i]
  r     = norm(buf.ri)
  E     = 0.5 * K * (r - req)^2
  @. buf.fi = -K * (r - req) * buf.ri / r

  F[i] .-= buf.fi
  F[j] .+= buf.fi

  E
end

function _harmonicBond!(
  F::Vf, u::Vu, lat::AbstractMatrix, i::Int64, j::Int64,
  K::Float64, req::Float64, rc::Float64;
  buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  r = pbcVec!(buf.ri, u[i], u[j], rc, lat)
  E = 0.5 * K * (r - req)^2

  @. buf.fi = -K * (r - req) * buf.ri / r

  F[i] .-= buf.fi
  F[j] .+= buf.fi

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
  F::Vf, u::Vu, i::Int64, o::Int64, j::Int64, K::Float64, θeq::Float64;
  buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  @. buf.ri = u[i] - u[o]
  @. buf.rj = u[j] - u[o]
  θ     = dot(buf.ri, buf.rj) / (norm(buf.ri) * norm(buf.rj)) |> (x -> clamp(x, -1, 1)) |> acos
  E     = 0.5 * K * (θ - θeq)^2
  pre   = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(buf.ri) * norm(buf.rj))
  buf.fi .= pre * (buf.rj - (buf.ri * (dot(buf.ri, buf.rj) / dot(buf.ri,buf.ri))))
  buf.fj .= pre * (buf.ri - (buf.rj * (dot(buf.ri, buf.rj) / dot(buf.rj,buf.rj))))
  F[i] .+= buf.fi
  F[j] .+= buf.fj
  @. F[o] -= (buf.fi + buf.fj)

  E
end

function _harmonicBondAngle!(
  F::Vf, u::Vu, lat::AbstractMatrix, i::Int64, o::Int64, j::Int64,
  K::Float64, θeq::Float64, rc::Float64;
  buf=_func_buffs
) where {Vf <: AbstractVector, Vu <: AbstractVector}

  di  = pbcVec!(buf.ri, u[o], u[i], rc, lat)
  dj  = pbcVec!(buf.rj, u[o], u[j], rc, lat)
  θ   = dot(buf.ri, buf.rj) / (di * dj) |> (x -> clamp(x, -1, 1)) |> acos
  E   = 0.5 * K * (θ - θeq)^2
  pre = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * di * dj)

  buf.fi .= pre * (buf.rj - (buf.ri * (dot(buf.ri, buf.rj) / dot(buf.ri,buf.ri))))
  buf.fj .= pre * (buf.ri - (buf.rj * (dot(buf.ri, buf.rj) / dot(buf.rj,buf.rj))))
  F[i]  .+= buf.fi
  F[j]  .+= buf.fj
  F[o]  .-= (buf.fi .+ buf.fj)

  E
end