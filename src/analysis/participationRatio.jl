
function getIPR(modes; N=2)
  m   = N * 3
  l   = size(modes)[2]
  n   = div(l, m)
  ipr = zeros(n,l)  

  for i = 1:l
    ρ = reshape(modes[:, i], (m, n))
    b = dot(modes[:, i], modes[:, i])^2
    
    for j = 1:n
      v         = ρ[:, j]
      ipr[j, i] = dot(v,v)^2 / b
    end

  end
  
  ipr
end

function getPR(ipr; x=1.6e-5)
  l  = size(ipr)[2]
  pr = zeros(l)

  for i = 1:l
    n     = findall(e -> e > x, ipr[:, i]) |> length
    pr[i] = n
  end

  pr
end