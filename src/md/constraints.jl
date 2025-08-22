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

struct SHAKE_BUF{F<:Float64, RO<:AbstractVector, AV3D<:AbstractVector}
  req::F
  rold::RO
  rnew::AV3D
end

SHAKE(bonds, req) = StaticConstraint(
  bonds, 
  SHAKE!,
  SHAKE_BUF(
    req,
    [zeros(3) for i = 1:length(bonds)],
    MVector{3}(zeros(3))
  ) 
)

function SHAKE!(
  F::AbstractVector, u::AbstractVector, m::AbstractVector,
  bonds::AbstractVector, buff::B
) where {B}

  maxIter = 25

  for (n,(i,j)) in enumerate(bonds)

    μ = 1 / ((1/m[i]) + (1/m[j]))
    
    # Set x to skip while loop if this is the first step
    x = if iszero(buff.rold[n])
      @. buff.rnew = u[j] - u[i]
      1e-9
    else
      1
    end

    iter = 0
    while abs(x) > 1e-8 && iter < maxIter

      @. buff.rnew = u[j] - u[i]

      x = 0.5 * (buff.req^2 - dot(buff.rnew, buff.rnew)) / dot(buff.rold[n], buff.rnew)

      u[i] .-= x * μ / m[i] .* buff.rold[n]
      u[j] .+= x * μ / m[j] .* buff.rold[n]

      iter += 1
    end

    buff.rold[n] .= buff.rnew
  end

end