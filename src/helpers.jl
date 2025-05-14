
#Calculate the shift frequency for a molecule at vib energy E
freqShiftMorse(v0, D, E) = v0 * ( (D - E) / D )^0.5

#Calculate the energy for a given shifted frequency
engyShiftMorse(v0, D, v) = D - D * (v/v0)^2

function mkvar(x)
  fn = Symbol(x)
  X  = @eval $fn
  return X
end

function myRepeat(A::Vector{T}, count::Integer, mask::Vector{Bool}) where T <: Union{Char, Number}
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

function getFrame(tj, i::Int64)
  m = tj.m
  s = tj.s
  r = tj.r[i]
  v = tj.v[i]

  [Atom(r[j], v[j], m[j], s[j]) for j in 1:length(m)]
end

function getFrame!(bdys, tj, i::Int64)
  r = tj.r[i]
  v = tj.v[i]

  for j = 1:length(bdys)
    bdys[j].r = r[j]
    bdys[j].v = v[j]
  end

end

function getLastFrame(solu)
  n = length(solu.prob.p.bdys)
  
  new = Atom[]
  for i in 1:n
    r = solu.u[end].x[2][i]
    v = solu.u[end].x[1][i]
    m = solu.prob.p.bdys[i].m
    s = solu.prob.p.bdys[i].s
    push!(new, Atom(r, v, m, s))
  end
  return new
end

function getLastFrame!(bdys::Vector{MyAtoms}, solu)
  N = length(bdys)
  
  for i in 1:N
    bdys[i].r = solu.u[end].x[2][i]
    bdys[i].v = solu.u[end].x[1][i]
    bdys[i].m = solu.prob.p.bdys[i].m
    bdys[i].s = solu.prob.p.bdys[i].s
  end
end

function getLastFrame!(cell::MyCell, solu)
  x0 = [j for i in solu.u[end].x[2] for j in i]

  cell.velocity   .= solu.u[end].x[1]
  cell.scaled_pos .= getScaledPos(x0, cell.lattice)

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
    bdys[i].v -= (x / bdys[i].m)
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

function getForces(EoM, bdys::Vector{MyAtoms})
  x0, vars = prep4pot(EoM, bdys)
  G        = zero(x0)
  EoM(nothing, G, x0, vars)

  [-G[i:i+2] for i = 1:3:length(G)]
end

function getForces(EoM, cell::MyCell)
  x0, vars = prep4pot(EoM, cell)
  G        = zero(x0)
  EoM(nothing, G, x0, vars)

  [-G[i:i+2] for i = 1:3:length(G)]
end