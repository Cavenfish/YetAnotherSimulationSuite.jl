"""
Constraints

Static Constraints should always have 
  apply!(forces, pos, mas, inds, buff)
  
Dynamic Constraints should always have
  apply!(forces, vels, pos, mas, inds buff)
"""

struct StaticConstraint{I,A,B} <: MyStaticConstraint
  inds::I
  apply!::A
  buff::B
end

fixAtoms(inds) = StaticConstraint(inds, fixAtoms!, nothing)

function fixAtoms!(forces, pos, mas, inds, buff)
  for i in inds
    forces[i] .= 0.0
  end
end

struct SHAKE_BUF{F<:Float64, AV3D<:AbstractVector}
  req::F
  rold::AV3D
  rnew::AV3D
end

SHAKE(bonds, req) = StaticConstraint(
  bonds, 
  SHAKE!,
  SHAKE_BUF(
    req,
    MVector{3}(zeros(3)),
    MVector{3}(zeros(3))
  ) 
)

function SHAKE!(
  F::AbstractVector, u::AbstractVector, m::AbstractVector,
  bonds::AbstractVector, buff::B
) where {B}

  maxIter = 25

  for (i,j) in bonds

    μ = 1 / ((1/m[i]) + (1/m[j]))    

    @. buff.rold = u[j] - u[i]

    x    = 1.0
    iter = 0
    while abs(x) > 1e-8 && iter < maxIter

      @. buff.rnew = u[j] - u[i]

      x = 0.5 * (buff.req^2 - dot(buff.rnew, buff.rnew)) / dot(buff.rold, buff.rnew)
      
      u[i] .-= x * μ / m[i] .* buff.rold
      u[j] .+= x * μ / m[j] .* buff.rold

      iter += 1
    end

    @. buff.rold = u[j] - u[i]
    x    = 1.0
    iter = 0
    while abs(x) > 1e-8 && iter < maxIter
      @. buff.rnew = F[j] / m[j] - F[i] / m[i]

      x = -dot(buff.rnew, buff.rold) / buff.req^2

      F[i] .-= x * μ .* buff.rold
      F[j] .+= x * μ .* buff.rold
    end
  end
end
