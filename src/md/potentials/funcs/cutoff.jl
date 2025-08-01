"""
Switching Functions for Potentials

"""

function switchSR!(S::AbstractVector, r::Float64, rs::Float64, rc::Float64)

  if r < rs
    S[1] = 1.0
    S[2] = 0.0
    return
  elseif r > rc
    S .= 0.0
    return
  end

  deno = (rc - rs)^3
  S[1] = (rc - r)^2 * (rc + 2r - 3rs) / deno
  S[2] = 6 * (rc - r) * (rs - r) / deno 
end

function switchLR!(S::AbstractVector, r::Float64, rs::Float64, rc::Float64)

  if r < rs
    S[1] = 1.0
    S[2] = 0.0
    return
  elseif r > rc
    S .= 0.0
    return
  end

  x    = (r - rs) / (rc - rs)
  S[1] = 1 - 10x^3 + 15x^4 - 6x^5
  S[2] = (-30x^2 + 60x^3 - 30x^4) / (rc-rs)
end