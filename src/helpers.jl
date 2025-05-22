
function mkvar(x)
  fn = Symbol(x)
  X  = @eval $fn
  return X
end

function myRepeat(A::Vector{T}, count::Integer, mask::Vector{Bool}) where T <: Union{String, Number}
  [A; repeat(A[.!mask], count-1)]
end

# Hate how this function is, but I couldn't find another solution
# to the A = Vector{Vector{Float64}} issue.
# I checked a rewrite of this where no middle man (B) was needed
# the allocations were unchanged, since this is cleaner to read
# I'm keeping it.
function myRepeat(A::Vector{Vector{Float64}}, count::Integer, mask::Vector{Bool})
  B = [A; repeat(A[.!mask], count-1)]
  C = zero(B)

  for (b,c) in zip(B,C)
    c .= b
  end

  C
end

function CoM(bdys)
  M = sum([i.m for i in bdys])
  r = sum([i.m*i.r for i in bdys])
  return r ./ M
end

function CoM(pos,mas)
  M = sum(mas)
  r = sum([mas[i]*pos[i] for i in 1:length(mas)])
  return r ./ M
end

function vCoM(bdys)
  M = sum([i.m for i in bdys])
  v = sum([i.m*i.v for i in bdys])
  return v ./ M
end

function vCoM(vel,mas)
  M = sum(mas)
  v = sum([mas[i]*vel[i] for i in 1:length(mas)])
  return v ./ M
end

function zeroVCoM!(bdys)
  N    = length(bdys)
  M    = sum([i.m for i in bdys])
  vcom = vCoM(bdys)
  x    = vcom / N * M
  
  for i in 1:N
    bdys[i].v .-= (x / bdys[i].m)
  end
end

function reducedMass(bdys)
  μ = sum([i.m^-1 for i in bdys])^-1
  return μ
end

function reducedMass(mas::Vector{Float64})
  μ = sum([m^-1 for m in mas])^-1
  return μ
end

function getForces(calc::MyCalc, bdys::Vector{MyAtoms})
  u       = [i.r for i in bdys]
  F       = zero(u)
  _, vars = prep4pot(calc.b, bdys)

  if calc.f! != nothing
    calc.f!(F, u, vars)
  else
    calc.ef!(F, u, vars)
  end

  F
end
  
function getForces(calc::MyCalc, cell::MyCell)
  u       = getPos(cell)
  F       = zero(u)
  _, vars = prep4pot(calc.b, cell)

  if calc.f! != nothing
    calc.f!(F, u, vars)
  else
    calc.ef!(F, u, vars)
  end

  F
end