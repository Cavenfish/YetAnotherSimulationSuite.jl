"""
Intramolecular Potential Functions
"""


function _Morse(r::Float64, rvec::Vector{Float64}, D::Float64, 
                a::Float64, req::Float64)
  c = exp(-a*(r-req))
  E = D * (1 - c)^2
  F = @. -2D * a * c * (1 - c) * rvec / r

  E, F
end

function _Morse!(F::Vector{Vector{Float64}}, u::Vector{Vector{Float64}},
                 i::Int64, j::Int64, D::Float64, a::Float64, req::Float64)
  
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  c     = exp(-a*(r-req))
  E     = D * (1 - c)^2
  f     = @. -2D * a * c * (1 - c) * rvec / r

  F[i] .-= f
  F[j] .+= f

  E
end

function _harmonicBond(r::Float64, rvec::Vector{Float64}, K::Float64, req::Float64)
  E     = 0.5 * K * (r - req)^2
  f     = @. - K * (r - req) * revc / r

  E, F
end

function _harmonicBond!(F::Vector{Vector{Float64}}, u::Vector{Vector{Float64}}, 
                        i::Int64, j::Int64, K::Float64, req::Float64)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  E     = 0.5 * K * (r - req)^2
  f     = @. - K * (r - req) * rvec / r

  F[i] .-= f
  F[j] .+= f

  E
end

function _harmonicBondAngle(r1::Vector{Float64}, r2::Vector{Float64}, 
                            K::Float64, θeq::Float64)
  θ   = dot(r1, r2) / (norm(r1) * norm(r2)) |> (x -> round(x, digits=10)) |> acos
  E   = 0.5 * K * (θ - θeq)^2
  pre = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(r1) * norm(r2))
  F1  = pre * (r2 - (r1 * (dot(r1, r2) / dot(r1,r1))))
  F2  = pre * (r1 - (r2 * (dot(r1, r2) / dot(r2,r2))))
  Fo  = @. - (F1 + F2)

  E, F1, F2, Fo
end

function _harmonicBondAngle!(F::Vector{Vector{Float64}}, u::Vector{Vector{Float64}}, 
                             i::Int64, o::Int64, j::Int64, K::Float64, θeq::Float64)
  ri    = u[i] - u[o]
  rj    = u[j] - u[o]
  θ     = dot(ri, rj) / (norm(ri) * norm(rj)) |> (x -> round(x, digits=10)) |> acos
  E     = 0.5 * K * (θ - θeq)^2
  pre   = K * (θ - θeq) / (sqrt(1 - cos(θ)^2) * norm(ri) * norm(rj))
  Fi    = pre * (rj - (ri * (dot(ri, rj) / dot(ri,ri))))
  Fj    = pre * (ri - (rj * (dot(ri, rj) / dot(rj,rj))))
  F[i] .+= Fi
  F[j] .+= Fj
  F[o] .-= (Fi + Fj)

  E
end