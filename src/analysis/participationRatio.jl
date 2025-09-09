"""
    getIPR(modes; N=2)

Compute the inverse participation ratio (IPR) for vibrational modes.

# Arguments
- `modes`: Matrix of mode vectors.
- `N`: Number of atoms per molecule (default: 2).

# Returns
- Matrix of IPR values.
"""
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

"""
    getPR(ipr; x=1.6e-5)

Compute the participation ratio (PR) from IPR data.

# Arguments
- `ipr`: IPR matrix.
- `x`: Threshold value (default: 1.6e-5).

# Returns
- Vector of participation ratios.
"""
function getPR(ipr; x=1.6e-5)
  l  = size(ipr)[2]
  pr = zeros(l)

  for i = 1:l
    n     = findall(e -> e > x, ipr[:, i]) |> length
    pr[i] = n
  end

  pr
end