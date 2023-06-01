
function vacf(vel, mas, dt; atms=false, norm=true)
  T = length(vel)       # Total timesteps 
  N = length(vel[1])    # Total number of atoms
  D = length(vel[1][1]) # Total number of dimensions
  c = zeros(Float64, T-1)
  
  for i in 1:N
    for j in 1:D
      data = [a[i][j] for a in vel]
      forw = fft(data)
      tmp  = forw .* conj(forw)
      back = ifft(tmp) ./ T
      c  .+= mas[i] * real(back[1:T-1])

  #Diffusion Coefficient
  

end
