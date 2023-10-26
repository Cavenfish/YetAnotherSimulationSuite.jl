
function Maxwell_Boltzmann(v, m, T, kB)
  #Need to make this for Langevin
end

function getTemp(m, v, kB, N)
  # Leave out the 1/2 to get 2Ekin for T calc
  Nf   = 3N - 3
  v2   = [i'i for i in v]
  Ekin = sum(m .* v2)
  Tsim = Ekin / (kB * Nf)
  return Tsim
end

struct Berendsen
  T::Float64
  kB::Float64
  gamma::Float64
end

function Berendsen!(T, a, v, m, inp)
  N    = length(m)
  Tsim = getTemp(m, v, inp.kB, N)
  push!(T, Tsim)

  if Tsim == 0.0
    a .+= inp.gamma .* v
  else
    a .+= inp.gamma * (inp.T / Tsim - 1) .* v
  end
end

struct Langevin
  T::Float64
  kB::Float64
  gamma::Float64
end

function Langevin!(T, a, v, m, inp)
  N    = length(m)
  Tsim = getTemp(m, v, inp.kB, N)
  push!(T, Tsim)

  # M   = sum(m) / N
  # sig = sqrt(2 * inp.gamma * M * inp.kB * inp.T)
  G   = Normal(0.0, 1)
  eta = [rand(G, 3) for i in 1:N] 

  sigma = @. sqrt(2 * inp.gamma * m * inp.kB * inp.T) 
  
  a1 = @. -inp.gamma * v
  a2 = @. sigma / m * eta

  @. a += a1 + a2
end

struct BDP
  T::Float64
  kB::Float64
  tau::Float64
end

struct BDPnT
  T::Vector
  mols::Vector
  kB::Float64
  tau::Float64
end

function BDPnT!(T, a, v, m, inp)
  for i in 1:length(inp.T)
    try 
      T[i] 
    catch err
      push!(T, [])
    end
    j   = inp.mols[i]
    tmp = BDP(inp.T[i], inp.kB, inp.tau)
    @views BDP!(T[i], a[j], v[j], m[j], tmp)
  end
end

function BDP!(T, a, v, m, inp)
  N    = length(m)
  Tsim = getTemp(m, v, inp.kB, N)
  push!(T, Tsim)

  if Tsim == 0.0
    return
  end

  Nf = (3N) - 3
  n  = Normal(0.0,1.0)
  R1 = rand(n)
  if (Nf-1)%2 == 0
    alpha = (Nf - 1)/2
    G     = Gamma(alpha, 1)
    R2    = 2 .* rand(G)
  else
    alpha = (Nf - 2)/2
    G     = Gamma(alpha, 1)
    R2    = 2 .* rand(G) .+ (rand(n) .^2)
  end

  if inp.tau > 0.1
    c1 = exp(-1/inp.tau)
  else
    c1 = 0.0
  end

  v2    = [i'i for i in v]
  K     = sum(0.5 .* m .* v2)
  sigma = 0.5 * Nf * inp.kB * inp.T

  Knew   = @. ( 
                K + (1-c1) * (sigma * (R2 + R1^2) / Nf - K) + 
                2 * R1 * sqrt(K * sigma / Nf * (1-c1) * c1)
              )
  alpha2 = @. sqrt(Knew / K)

  @. v *= alpha2
end