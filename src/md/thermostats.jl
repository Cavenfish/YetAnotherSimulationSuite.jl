
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

function Berendsen!(T, a, v, m, inp)
  N    = length(m)
  tmp  = (m/2) .* v 
  Ekin = tmp'tmp
  Tsim = Ekin / (inp.kB * (3 * N))
  push!(T, Tsim)

  if Ekin == 0.0
    a .+= inp.gamma .* v
  else
    a  .+= inp.gamma * (inp.T / Tsim - 1) .* v
  end
end

function Langevin!(T, a, v, m, inp)
  tmp  = (m/2) .* v 
  Ekin = tmp'tmp
  Tsim = Ekin / (inp.kB * (3 * N))
  push!(T, Tsim)

  a1 = -inp.gamma .* v ./ m

  