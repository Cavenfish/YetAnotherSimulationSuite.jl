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
"""

function g(i,j)
  ita  = 0.3
  diff = j - i
  r    = sqrt(diff'diff)
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
  @views P[3] = sqrt(P[3])
  @views P[4] = sqrt(P[4])
  @views P[5] = sqrt(P[5])

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

function readInVars(file)
  """
  What a nightmare it is to read this in. 

  Moving forward I will get rid of this function. 
  Instead I will store vars in faster reading format.
  """
  s0 = 7
  s1 = 45
  s2 = 45
  w1 = zeros(Float64, s0, s1)
  b1 = zeros(Float64, s1)
  w2 = zeros(Float64, s1, s2)
  b2 = zeros(Float64, s2)
  w3 = zeros(Float64, s2)
  b3 = 0.0
  rg = zeros(Float64, 2, s0)
  vg = zeros(Float64, 2)


  inp = readlines(file)

  # The txt file is split by blank lines
  # where each block of numbers (number between blanks)
  # are for one variable. 
  x = Int32[]
  push!(x, 1)
  for i in findall(e -> e  == " ", inp)
    push!(x, i-1)
    push!(x, i+1)
  end

  # Here we are filling in the vars similar to how
  # the original code does it. We fill each var in
  # based on the numbers in its block. We fill in 
  # column-wise (ie. w1[1st, 2nd]).

  w1[:,:] = parse.(Float64, inp[x[1]  : x[2]])
  b1[:]   = parse.(Float64, inp[x[3]  : x[4]])
  w2[:,:] = parse.(Float64, inp[x[5]  : x[6]])
  b2[:]   = parse.(Float64, inp[x[7]  : x[8]])
  w3[:]   = parse.(Float64, inp[x[9]  : x[10]])
  b3      = parse.(Float64, inp[x[11] : x[12]])[1] #[1] to unvector it

  # Only these last two block are double column 
  # why? because making it easy would be dumb
  tmp = inp[x[13] : x[14]]
  tmp = split.(tmp, " ")
  deleteat!.(tmp, findall.(e -> e == "", tmp))

  for i in 1:s0
    rg[:,i] = parse.(Float64, tmp[i])
  end

  tmp   = inp[x[15] : x[16]][1]
  tmp   = split(tmp, " ")
  deleteat!(tmp, findall(e -> e == "", tmp))
  vg[:] = parse.(Float64, tmp)

  vgg = vg[2] - vg[1]
  rgg = rg[2,:] - rg[1,:]

  return w1,b1,w2,b2,w3,b3,rg,rgg,vg,vgg
end

function pairPot(co1, co2, vars, dPdr, P)
  # Get PIPs: P is a vector 
  # P, dPdr = getPIPs(co1..., co2...)
  getPIPs!(P,dPdr,co1..., co2...)

  # Weights and biases 
  w1,b1,w2,b2,w3,b3,rg,rgg,vg,vgg = vars

  # Map min-max
  @views P[:] = 2 * (P[:]-rg[1,:]) ./ rgg[:] .- 1

  # 1st layer
  y = b1 + transpose(w1) * P
  f = tanh.(y)
  A = w1 * (w2 .* (1 .- f.^2))

  # 2nd layer
  y = b2 + transpose(w2) * f
  f = tanh.(y)

  # output layer
  v  = b3 + w3'f
  dv = A * (w3 .* (1 .- f.^2))

  # remapping
  v   = vgg * (v+1)/2 + vg[1]
  dv  = (dv .* vgg) ./ rgg

  dv  = transpose(dPdr) * dv
  return v, dv
end

function molPot(mol)
  # Define weights and biases
  w1 = [2.5731334302226245e0,
       -4.1635475296365723e0,
        6.2652688635816842e0,
       -8.1228824615517947e0,
        2.5634824563114837e+1,
       -1.6643878666848999e0,
        1.2863481593265147e+1,
       -5.1395051186435685e0]

  b1 = [-1.6421783363943476e0,
         1.5364774081764951e0,
        -1.1338455332512607e0,
        -1.4448051810696709e0,
         7.5217573991947644e0,
        -1.4005229119446290e0,
         1.1053854378210930e+1,
        -5.9299626180269485e0]

  w2 = [-8.8664820626030757e-3,
         7.8571245473773067e-3,
        -3.6411047563733342e-3,
        -4.0358215533209145e-3,
         9.6640587889626451e-4,
        -1.4325782866595651e0,
         1.2002568907875554e-2,
         8.3983298757280007e0]
  
  b2 = 6.8970075393140338
  ra = 1.4
  rb = 7.0
  va = 9.8087308941410573e-2
  vb = 1.9558422718340193e+5

  # Now we can start the Pot Calc
  
  # Get bond length (r)
  diff = mol[2] - mol[1]
  r    = sqrt(diff'diff)
  rhat = diff / r

  # Map min-max
  x  = 2.0 * (r-ra)/(rb-ra) - 1.0

  # 1st layer
  y  = b1 + w1 * x
  f  = tanh.(y)

  # output layer
  v  = b2 + w2'f
  dv = (1 .- f.^2) .* w1
  dv = w2'dv

  # remapping
  v   = (v+1) * (vb-va)/2.0 + va
  v  += 0.560096850315234
  dv *= (vb-va) / (rb-ra)

  # apply rhat to get F
  F = dv * rhat
  return v, dv, F
end

function getUnitVectors!(r, co1, co2)
  rhat(v) = v / sqrt(v'v)

  c1,o1 = co1
  c2,o2 = co2

  r[1,:] = rhat(c1 - o1)
  r[2,:] = rhat(c2 - o2)
  r[3,:] = rhat(c2 - o1)
  r[4,:] = rhat(o2 - c1)
  r[5,:] = rhat(o1 - o2)
  r[6,:] = rhat(c1 - c2)
end

function HGNNdyn(a, v, u, p, t)

  # initialize things
  E = 0.0
  F = zero(u)
  m = [i.m for i in p.bdys]
  r = u ./ 0.5291772083 # to Bohr
  
  # Pre-allocate for performance gains
  rhats = zeros(Float64, 6, 3)
  dPdr  = zeros(Float64, 7, 6)
  P     = zeros(Float64, 7)
  
  # Get weight and biases:
  #   - weights are matricies (except w3)
  #   - biases are vectors (except b3)
  # For now, I will hardcode the input file
  inp  = "/home/brian/Research/JMD/ogSRC/nn_ococ_w20.txt"
  vars = readInVars(inp)

  for i in p.mols
    v, dv, f = molPot(r[i])
    E       += v
    F[i[1]] += f
    F[i[2]] -= f
  end

  for i in p.pars
    c1,o1  = i[1]
    c2,o2  = i[2]
    v, dv  = pairPot(r[i[1]], r[i[2]], vars, dPdr, P)
    E     += v
    @views getUnitVectors!(rhats, r[i[1]], r[i[2]])
    rhats .*= dv

    # rhat: o1 --> c1
    @views F[o1] += rhats[1,:]
    @views F[c1] -= rhats[1,:]
    
    # rhat: o2 --> c2
    @views F[o2] += rhats[2,:]
    @views F[c2] -= rhats[2,:]

    # rhat: o1 --> c2
    @views F[o1] += rhats[3,:]
    @views F[c2] -= rhats[3,:]

    # rhat: c1 --> o2
    @views F[c1] += rhats[4,:]
    @views F[o2] -= rhats[4,:]

    # rhat: o2 --> o1
    @views F[o2] += rhats[5,:]
    @views F[o1] -= rhats[5,:]

    # rhat: c2 --> c1
    @views F[c2] += rhats[6,:]
    @views F[c1] -= rhats[6,:]
  end
  
  E  *= 0.000124 # cm-1 to eV
  F .*= (0.000124 / 0.5291772083) # cm-1/Bohr to eV/Angstrom

  a .= F ./ m
  push!(p.time, t)
  push!(p.energy, E)
  push!(p.forces, F)
end

function HGNNpot(F, G, y0, p)

  # Pre-allocate for performance gains
  rhats = zeros(Float64, 6, 3)
  dPdr  = zeros(Float64, 7, 6)
  P     = zeros(Float64, 7)

  # I couldn't get Optim to work with a 2D vector
  # so I had to flatten the vector before sending 
  # it through. This worked but I then need to 
  # flatten -> send -> unflatten. Which is such 
  # a pain. I need a better solution.
  x0     = Vector[]
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(x0, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

  # initialize things
  energy = 0.0
  r = x0 ./ 0.5291772083 # to Bohr
  
  # Get weight and biases:
  #   - weights are matricies (except w3)
  #   - biases are vectors (except b3)
  # For now, I will hardcode the input file
  inp  = "/home/brian/Research/JMD/ogSRC/nn_ococ_w20.txt"
  vars = readInVars(inp)

  for i in p.mols
    v, dv, f = molPot(r[i])
    energy  += v
    forces[i[1]] += f
    forces[i[2]] -= f
  end

  for i in p.pars
    c1,o1   = i[1]
    c2,o2   = i[2]
    v, dv   = pairPot(r[i[1]], r[i[2]], vars, dPdr, P)
    energy += v
    @views getUnitVectors!(rhats, r[i[1]], r[i[2]])
    rhats .*= dv

    # rhat: o1 --> c1
    @views forces[o1] += rhats[1,:]
    @views forces[c1] -= rhats[1,:]
    
    # rhat: o2 --> c2
    @views forces[o2] += rhats[2,:]
    @views forces[c2] -= rhats[2,:]

    # rhat: o1 --> c2
    @views forces[o1] += rhats[3,:]
    @views forces[c2] -= rhats[3,:]

    # rhat: c1 --> o2
    @views forces[c1] += rhats[4,:]
    @views forces[o2] -= rhats[4,:]

    # rhat: o2 --> o1
    @views forces[o2] += rhats[5,:]
    @views forces[o1] -= rhats[5,:]

    # rhat: c2 --> c1
    @views forces[c2] += rhats[6,:]
    @views forces[c1] -= rhats[6,:]
  end
  
  energy  *= 0.000124 # cm-1 to eV
  forces .*= (0.000124 / 0.5291772083) # cm-1/Bohr to eV/Angstrom

  if G != nothing
    tmp = [j for i in forces for j in i]
    for i in 1:length(G)
      G[i] = -tmp[i]
    end
  end

  if F != nothing
    return energy
  end
end