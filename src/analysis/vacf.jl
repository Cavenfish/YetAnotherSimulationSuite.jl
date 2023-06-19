
struct vacfInps
  vel::Vector
  mas::Vector{Float64}
  Hz::Float64
  norm::Bool
  win::Function
  pad::UInt8
  mir::Bool
end

mutable struct vacfOut
  c::Vector{Float64}
  C::Vector{Float64}
  D::Float64
  v::Vector{Float64}
  I::Vector{Float64}
end

Hann( i,N) = sin((pi*i)/N)^2
Welch(i,N) = 1 - ((i-N/2)/(N/2))^2
HannM(i,N) = cos((pi*i)/(2*(N-1)))^2

function window!(out, W)
  N = length(out.C)
  w = [W(i,N) for i in 0:N-1]

  out.C .*= w
end

function padZeros!(out, f)
  N = length(out.C)
  z = zeros(N*f)

  @views z[1:N] = out.C
  out.C = z
end

function mirror!(out)
  out.C = [reverse(out.C); out.C[2:end]]
end

function vacf!(inp, out; atms=false)
  T = length(inp.vel)       # Total timesteps 
  N = length(inp.vel[1])    # Total number of atoms
  D = length(inp.vel[1][1]) # Total number of dimensions
  
  out.c = zeros(Float64, T-1)
  
  for i in 1:N
    for j in 1:D
      data = [a[i][j] for a in inp.vel]
      forw = fft(data)
      tmp  = forw .* conj(forw)
      back = ifft(tmp) ./ T
      
      @views out.c  .+= inp.mas[i] * real(back[1:T-1])
    end
  end

  # out.D  = sum(c) / (D*N)
  # out.D *= UNIT  CONVERSION
  
  if inp.norm
    out.C = out.c ./ out.c[1]
  else
    out.C = out.c
  end

end

function VDOS(inp)
  dum = zeros(10)
  out = vacfOut(dum,dum,0.0,dum,dum)

  vacf!(inp, out)
  window!(out, inp.win)
  padZeros!(out, inp.pad)
  
  if inp.mir
    mirror!(out)
  end

  k = 29979245800.0 # Divide by this to convert from Hz to cm^-1
  N = length(out.C)
  I = fft(out.C)
  v = fftfreq(N, inp.Hz) ./ k
  
  n::UInt32 = div(length(v), 2)
  @views I  = I[1:n]
  @views v  = v[1:n]

  out.v = abs.(v)
  out.I = abs.(I)

  return out
end

function getVelMas(solu)
  T   = length(solu.t)
  vel = [solu.u[i].x[1] for i in 1:T]
  mas = [i.m for i in solu.prob.p.bdys]
  
  return vel, mas
end
