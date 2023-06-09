
function Maxwell_Boltzmann(v, m, T, kB)
  #Need to make this for Langevin
end

function getTemp(m, v, kB, N)
  tmp  = (m) .* v # Leave out the 1/2 to get 2Ekin for T calc
  Ekin = tmp'tmp
  Tsim = Ekin / (kB * (3 * N))
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
  dt::Float64
  tau::Float64
end

function BDP!(T, a, v, m, inp)
  N    = length(m)
  Tsim = getTemp(m, v, inp.kB, N)
  push!(T, Tsim)

  if Tsim == 0.0
    return
  end

  n = Normal(0.0,1.0)
  if (3N-1)%2 == 0
    alpha = (3N - 1)/2
    G     = Gamma(alpha, 1)
  else
    alpha = (3N - 2)/2
    G     = Gamma(alpha, 1)
  end

  R1 = rand(n, N) .^2
  R2 = rand(G, N)

  tmp  = (m/2) .* v 
  K    = (tmp'tmp) / N
  Kbar = 3N * inp.kB * inp.T

  c1 = exp(- inp.dt / inp.tau)
  c2 = Kbar / (3N * K)
  c3 = exp(- inp.dt / (2*inp.tau))

  alpha2 = @. c1 + c2 * (1-c1) * (R1 + R2) + 2*c3 * sqrt(c2 * (1-c1) * R1)

  @. v *= sqrt(alpha2) / m
end