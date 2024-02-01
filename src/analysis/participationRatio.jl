
function getIPR(modes; N=2)
  m   = N * 3
  l   = size(modes)[2]
  ipr = zeros(l)

  for i = 1:l
    ρ = reshape(modes[:, i], (m, l))
    a = 0.0
    b = 0.0

    for j = 1:l
      v  = ρ[:, j]
      a += dot(v, v)^2
      b += dot(v, v)
    end
    
    ipr[i] += a / (b^2)
  end
  
  ipr
end