"""
CO-CO Potential from Chen 2020

Notation
  - Symmetry functions Vector => P
  - Transforming equations    => g
  - Weights                   => w
  - Bias                      => b

Units
  - Energy  --> cm-1
  - Spatial --> Bohr

Author: Brian C. Ferrari

------ My Notes While Coding ------

dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1) 
    -> 1.d0 * Transpose(w1) * r + 1.d0 * rt1

dgemm -> matrix -- matrix multiplication 

dpdr takes the derivative of P wrt g. It becomes
a matrix because it does dpdr[:,1] = del/delr1 P[:],
where r1 =  r_{o1, c1}. It then follows like that
for the rest. This can be seen as simple vector math
Essentially:
    dPdr = P[7] * dr[6]
     ^      ^       ^
  matrix  column   row
          vector  vector 

Switching to PotVars struct has changed allocations in two ways
  1) Significantly increased allocations during for calls
     to the potential. This comes from loading the molVars
     and pairVars during the call rather than at module load.
  2) Significantly decreased allocations associated with the
     following vars: A, P, dPdr, rhats. This comes from them
     acting like a cache now rather than being reallocated at
     each new call.

The overall result is that for quick uses of the potential
(ie. single point energy) this method is slower and allocates
more than the previous method. However, for longer uses of the
potential (ie. md simulations) this significantly speeds up the
code by reducing allocations.
"""

HGNN() = Calculator(HGNN; EF=HGNN!)

struct _HGNN_PotVars{MV,PV,D1,D2,D3, F<:AbstractFloat} <: PotVars
  A::Matrix{F}
  P::SizedVector{D1,F}
  dPdr::SizedMatrix{D1,D2,F}
  rhats::SizedMatrix{D2,D3,F}
  molVars::MV
  pairVars::PV
end 

function HGNN(bdys::Vector{MyAtoms})
  f = jldopen(joinpath(@__DIR__, "params/hgnn.jld2"))
  molVars  = f["molVars"]
  pairVars = f["pairVars"]
  close(f)

  # Pre-allocate for performance gains
  rhats = SizedMatrix{6,3}(zeros(Float64, 6, 3))
  dPdr  = SizedMatrix{7,6}(zeros(Float64, 7, 6))
  P     = SizedVector{7}(zeros(Float64, 7))
  A     = zeros(Float64, 7, 45)

  _HGNN_PotVars(A, P, dPdr, rhats, molVars, pairVars)
end

function HGNN!(F, u, p)
  E = 0.0
  P = p.potVars
  r = u ./ 0.5291772083 # to Bohr

  for i in p.mols
    E += molPot!(F, r, i, P.molVars)
  end

  for i in p.pars
    E += pairPot!(F, r, i, P.pairVars, P.rhats, P.dPdr, P.P, P.A)
  end
  
  E  *= 0.000124 # cm-1 to eV
  F .*= (0.000124 / 0.5291772083) # cm-1/Bohr to eV/Angstrom

  E
end

function g(i,j)
  ita  = 0.3
  diff = j - i
  r    = sqrt(dot(diff,diff))
  p    = exp(-r * ita)
  return p
end

function getPIPs!(P,dPdr,c1,o1,c2,o2)

  # Predefine possible g calls to prevent repeating math
  g1 = g(o1, c1)
  g2 = g(o2, c2)
  g3 = g(o1, c2)
  g4 = g(o2, c1)
  g5 = g(o1, o2)
  g6 = g(c1, c2)

  # See Chen et al. 2020 for written out PIPs
  P[1] = g1   + g2
  P[2] = g3   + g4
  P[3] = g1^2 + g2^2
  P[4] = g3^2 + g4^2
  P[5] = g1 * g3 + g4 * g2
  P[6] = g5
  P[7] = g6

  #This part is in their code but not the paper
  P[3] = sqrt(P[3])
  P[4] = sqrt(P[4])
  P[5] = sqrt(P[5])

  # Note: ita = 0.3
  dPdr[1,1] = -0.3*g1 # dp1/dr1
  dPdr[1,2] = -0.3*g2 # dp1/dr2
  dPdr[2,3] = -0.3*g3 # dp2/dr3
  dPdr[2,4] = -0.3*g4 # dp2/dr4
  dPdr[6,5] = -0.3*g5 # dp6/dr5
  dPdr[7,6] = -0.3*g6 # dp7/dr6

  dPdr[3,1] = (-0.3 * g1^2) / P[3] # dp3/dr1
  dPdr[3,2] = (-0.3 * g2^2) / P[3] # dp3/dr2 
  dPdr[4,3] = (-0.3 * g3^2) / P[4] # dp4/dr3 
  dPdr[4,4] = (-0.3 * g4^2) / P[4] # dp4/dr4 

  dPdr[5,1] = (-0.3*g1*g3) / (2*P[5]) # dp5/dr1
  dPdr[5,2] = (-0.3*g4*g2) / (2*P[5]) # dp5/dr1
  dPdr[5,3] = (-0.3*g1*g3) / (2*P[5]) # dp5/dr1
  dPdr[5,4] = (-0.3*g4*g2) / (2*P[5]) # dp5/dr1
end

function pairPot!(F, u, i, vars, rhats, dPdr, P, A)
  # ref i
  c1, o1, c2, o2 = i[1][1], i[1][2], i[2][1], i[2][2]
  
  # Get PIPs: P is a vector, dPdr is a matrix
  getPIPs!(P,dPdr, u[c1], u[o1], u[c2], u[o2])
  getUnitVectors!(rhats, u[c1], u[o1], u[c2], u[o2])

  # Weights and biases 
  w1,b1,w2,b2,w3,b3,rg,rgg,vg,vgg = vars

  # Map min-max
  @views P .= 2 * (P[:] - rg[1,:]) ./ rgg[:] .- 1

  # 1st layer
  y = b1 + transpose(w1) * P
  f = tanh.(y)

  #Zero out matrix, then inplace fill with multiplied matrices 
  A .= 0.0
  mul!(A, w1, @. (w2 * (1 - f^2)))

  # 2nd layer
  y .= b2 + transpose(w2) * f
  f .= tanh.(y)

  # output layer
  v  = b3 + dot(w3,f)
  dv = A * @. (w3 * (1 - f^2))

  # remapping
  v      = vgg * (v+1)/2 + vg[1]
  @. dv  = (dv * vgg) / rgg

  # make forces in rhat directions
  rhats .*= transpose(dPdr) * dv

  # rhat: o1 --> c1
  F[o1] .+= rhats[1,:]
  F[c1] .-= rhats[1,:]
  
  # rhat: o2 --> c2
  F[o2] .+= rhats[2,:]
  F[c2] .-= rhats[2,:]

  # rhat: o1 --> c2
  F[o1] .+= rhats[3,:]
  F[c2] .-= rhats[3,:]

  # rhat: c1 --> o2
  F[c1] .+= rhats[4,:]
  F[o2] .-= rhats[4,:]

  # rhat: o2 --> o1
  F[o2] .+= rhats[5,:]
  F[o1] .-= rhats[5,:]

  # rhat: c2 --> c1
  F[c2] .+= rhats[6,:]
  F[c1] .-= rhats[6,:]

  v
end

function molPot!(F, u, i, vars)
  # Define weights and biases
  w1,b1,w2,b2,ra,rb,va,vb = vars
  
  # # Get bond length (r)
  diff = u[i[2]] - u[i[1]]
  r    = sqrt(dot(diff,diff))
  rhat = diff / r

  # Map min-max
  x  = 2.0 * (r-ra)/(rb-ra) - 1.0

  # 1st layer
  y  = b1 + w1 * x
  f  = tanh.(y)

  # output layer
  v  = b2 + dot(w2,f)
  dv = dot(w2, @. (1 - f^2) * w1)

  # remapping
  v   = (v+1) * (vb-va)/2.0 + va
  v  += 0.560096850315234
  dv *= (vb-va) / (rb-ra)

  # Inplace update forces
  @. F[i[1]] += dv * rhat
  @. F[i[2]] -= dv * rhat

  v
end

function getUnitVectors!(r, c1, o1, c2, o2)
  rhat(v) = v / sqrt(dot(v,v))

  r[1,:] = rhat(c1 - o1)
  r[2,:] = rhat(c2 - o2)
  r[3,:] = rhat(c2 - o1)
  r[4,:] = rhat(o2 - c1)
  r[5,:] = rhat(o1 - o2)
  r[6,:] = rhat(c1 - c2)
end