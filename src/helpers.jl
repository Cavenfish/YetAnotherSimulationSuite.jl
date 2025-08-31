"""
    mkvar(x)

Evaluate a symbol from a string and return its value.

# Arguments
- `x`: String representing the variable name.

# Returns
- Value of the variable with name `x`.
"""
function mkvar(x)
  fn = Symbol(x)
  X  = @eval $fn
  return X
end

"""
    myRepeat(A::Vector{T}, count::Integer, mask::Vector{Bool}) where T <: Union{String, Number}

Repeat elements of a vector, excluding masked elements, and concatenate with the original.

# Arguments
- `A`: Input vector.
- `count`: Number of times to repeat unmasked elements.
- `mask`: Boolean mask vector (true = keep, false = repeat).

# Returns
- Concatenated vector with repeated elements.
"""
function myRepeat(A::Vector{T}, count::Integer, mask::Vector{Bool}) where T <: Union{String, Number}
  [A; repeat(A[.!mask], count-1)]
end

# Hate how this function is, but I couldn't find another solution
# to the A = Vector{Vector{Float64}} issue.
# I checked a rewrite of this where no middle man (B) was needed
# the allocations were unchanged, since this is cleaner to read
# I'm keeping it.
"""
    myRepeat(A::Vector{MVector{D,Float64}}, count::Integer, mask::Vector{Bool}) where D

Repeat elements of a vector of static vectors, excluding masked elements, and return a deep copy.

# Arguments
- `A`: Vector of `MVector{D,Float64}`.
- `count`: Number of times to repeat unmasked elements.
- `mask`: Boolean mask vector (true = keep, false = repeat).

# Returns
- Vector of repeated and copied static vectors.
"""
function myRepeat(A::Vector{MVector{D,Float64}}, count::Integer, mask::Vector{Bool}) where D
  B = [A; repeat(A[.!mask], count-1)]
  C = zero(B)

  for (b,c) in zip(B,C)
    c .= b
  end

  C
end

"""
    CoM(bdys)

Compute the center of mass for a collection of atoms.

# Arguments
- `bdys`: Vector of objects with fields `m` (mass) and `r` (position).

# Returns
- Center of mass as a 3-element vector.
"""
function CoM(bdys)
  o  = zeros(3)
  M  = sum([i.m for i in bdys])
  r  = sum([i.m*i.r for i in bdys])
  o .= r ./ M

  o
end

"""
    CoM(pos, mas)

Compute the center of mass from positions and masses.

# Arguments
- `pos`: Vector of position vectors.
- `mas`: Vector of masses.

# Returns
- Center of mass as a vector.
"""
function CoM(pos,mas)
  M = sum(mas)
  r = sum([mas[i]*pos[i] for i in 1:length(mas)])
  return r ./ M
end

"""
    vCoM(bdys)

Compute the center of mass velocity for a collection of atoms.

# Arguments
- `bdys`: Vector of objects with fields `m` (mass) and `v` (velocity).

# Returns
- Center of mass velocity as a vector.
"""
function vCoM(bdys)
  M = sum([i.m for i in bdys])
  v = sum([i.m*i.v for i in bdys])
  return v ./ M
end

"""
    vCoM(vel, mas)

Compute the center of mass velocity from velocities and masses.

# Arguments
- `vel`: Vector of velocity vectors.
- `mas`: Vector of masses.

# Returns
- Center of mass velocity as a vector.
"""
function vCoM(vel,mas)
  M = sum(mas)
  v = sum([mas[i]*vel[i] for i in 1:length(mas)])
  return v ./ M
end

"""
    zeroVCoM!(bdys)

Remove the center of mass velocity from a collection of atoms in-place.

# Arguments
- `bdys`: Vector of objects with fields `m` (mass) and `v` (velocity).

# Side Effects
- Modifies velocities in-place to set total momentum to zero.
"""
function zeroVCoM!(bdys)
  N    = length(bdys)
  M    = sum([i.m for i in bdys])
  vcom = vCoM(bdys)
  x    = vcom / N * M
  
  for i in 1:N
    bdys[i].v .-= (x / bdys[i].m)
  end
end

"""
    reducedMass(bdys)

Compute the reduced mass for a collection of atoms.

# Arguments
- `bdys`: Vector of objects with field `m` (mass).

# Returns
- Reduced mass (Float64).
"""
function reducedMass(bdys)
  μ = sum([i.m^-1 for i in bdys])^-1
  return μ
end

"""
    reducedMass(mas::Vector{Float64})

Compute the reduced mass from a vector of masses.

# Arguments
- `mas`: Vector of masses.

# Returns
- Reduced mass (Float64).
"""
function reducedMass(mas::Vector{Float64})
  μ = sum([m^-1 for m in mas])^-1
  return μ
end

"""
    getForces(calc::MyCalc, bdys::Vector{MyAtoms})

Compute the forces on a set of atoms using a calculator.

# Arguments
- `calc`: Calculator object.
- `bdys`: Vector of `MyAtoms` objects.

# Returns
- Vector of force vectors.
"""
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

"""
    getForces(calc::MyCalc, cell::MyCell)

Compute the forces on all atoms in a cell using a calculator.

# Arguments
- `calc`: Calculator object.
- `cell`: `MyCell` object.

# Returns
- Vector of force vectors.
"""
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