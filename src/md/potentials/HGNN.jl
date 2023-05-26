"""
CO-CO Potential from Chen 2020

Notation
  - Symmetry functions Vector => P
  - Transforming equations    => g
  - Weights                   => w
  - Bias                      => b

Author: Brian C. Ferrari
"""

"""
Notes to keeps in mind:

dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1) 
    -> 1.d0 * Transpose(w1) * r + 1.d0 * rt1
"""

function g(i,j)
  lambda = 0.3#???
  diff   = j - i
  r      = sqrt(diff'diff)
  p      = exp(-r * lambda)
  return p
end

function getPIPs(c1,o1,c2,o2)

  p1 = g(o1, c1)   + g(o2, c2)
  p2 = g(o1, c2)   + g(o2, c1)
  p3 = g(o1, c1)^2 + g(o2, c2)^2
  p4 = g(o1, c2)^2 + g(o2, c1)^2

  p5 = g(o1, c1) * g(o1, c2) + g(o2,  c1) * g(o2, c2)
  p6 = g(o1, o2)
  p7 = g(c1, c2)

  #This part is in their code but not the paper
  p3 = sqrt(p3)
  p4 = sqrt(p4)
  p5 = sqrt(p5)
  
  return [p1, p2, p3, p4, p5, p6, p7]
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

  return w1,b1,w2,b2,w3,b3,rg,vg
end

function pairPot(pair)
  # For now, I will hardcode the input file
  inp = "/home/brian/Research/JMD/ogSRC/nn_ococ_w20.txt"

  # Get PIPs: P is a vector 
  P = getPIPs(pair...)

  # Get weight and biases:
  #   - weights are matricies (except w3)
  #   - biases are vectors (except b3)
  w1,b1,w2,b2,w3,b3,rg,vg = readInVars(inp)

  # Map min-max
  rgg  = rg[2,:] - rg[1,:]
  P[:] = 2 * (P[:]-rg[1,:]) ./ rgg[:] .- 1

  # 1st layer
  y = b1 + transpose(w1) * P
  f = tanh.(y)

  # 2nd layer
  y = b2 + transpose(w2) * f
  f = tanh.(y)

  # output layer
  v = b3 + w3'f

  # remapping
  vgg = vg[2] - vg[1]
  v   = vgg * (v+1)/2 + vg[1]
  return v
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

  # Map min-max
  x = 2.0 * (r-ra)/(rb-ra) - 1.0

  # 1st layer
  y = b1 + w1 * x
  f = tanh.(y)

  # output layer
  v = b2 + w2'f

  # remapping
  v  = (v+1) * (vb-va)/2.0 + va
  v += 0.560096850315234
  return v
end


