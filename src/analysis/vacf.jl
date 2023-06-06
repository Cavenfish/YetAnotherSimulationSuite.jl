
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

function window!(c, W)
  N = length(c)

  for i in 1:N
    c[i] *= W(i,N)
  end
end

function padZeros!(c, f)
  N = length(c)
  z = zeros(N*f)

  @views z[1:N] = c
  
  c = z
end

function mirror!(c)
  c = [reverse(c); c[2:end]]
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
  end

end

function VDOS(inp)
  dum = zeros(10)
  out = vacfOut(dum,dum,0.0,dum,dum)

  vacf!(inp, out)
  window!(out.C, inp.win)
  padZeros!(out.C, inp.pad)
  
  if inp.mir
    mirror!(out.C)
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
