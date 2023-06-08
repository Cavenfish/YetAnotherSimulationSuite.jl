
struct Berendsen
  T::Float64
  kB::Float64
  gamma::Float64
end

struct Langevin
  T::Float64
  kB::Float64
  gamma::Float64
end

function Maxwell_Boltzmann(v, m, T, kB)
  #Need to make this for Langevin
end


function Berendsen!(T, a, v, m, inp)
  N    = length(m)
  tmp  = (m/2) .* v 
  Ekin = tmp'tmp
  Tsim = Ekin / (inp.kB * (3 * N))
  push!(T, Tsim)

  if Ekin == 0.0
    a .+= inp.gamma .* v
  else
    a .+= inp.gamma * (inp.T / Tsim - 1) .* v
  end
end

function Langevin!(T, a, v, m, inp)
  N    = length(m)
  tmp  = (m/2) .* v 
  Ekin = tmp'tmp
  Tsim = Ekin / (inp.kB * (3 * N))
  push!(T, Tsim)

  M   = sum(m) / N
  mu  = sqrt(3 * inp.kB * inp.T / M)
  sig = sqrt(inp.kB * inp.T / M)
  G   = Normal(mu, sig)
  eta = [rand(G, 3) for i in 1:N]

  sigma = @. sqrt(2 * inp.gamma * m * inp.kB * inp.T) 
  
  a1 = @. -inp.gamma * v
  a2 = @. sigma / m * eta
  

  @. a += a1 + a2
end


  