"""
Damping Functions for Potentials

Notes:
  - Functions will be passed to inter functions,
    so they should follow the format foo(r, p).
    Where p can be a single value or an array with
    all additional damping parameters.
"""

function tangToennies(r, D)
  dr = D * r
  s  = [dr^i/factorial(i) for i = 0:6] |> sum
  
  1 - (exp(-dr) * s)
end