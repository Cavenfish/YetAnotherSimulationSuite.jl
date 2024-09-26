"""
Damping Functions for Potentials

Notes:
  - Functions will be passed to inter functions,
    so they should follow the format foo(r, p).
    Where p can be a single value or an array with
    all additional damping parameters.
  - All functions must return the energy and force
    damping as a tuple like (f,fp).
"""

function tangToennies(r, D)
  dr = D * r
  s  = [dr^i/factorial(i) for i = 0:6] |> sum
  sp = [i*(dr)^(i-1)/factorial(i) for i = 0:6] |> sum
  
  f  = 1 - (exp(-dr) * s)
  fp = D * exp(-dr) * s - exp(-dr) * sp

  (f,fp)
end