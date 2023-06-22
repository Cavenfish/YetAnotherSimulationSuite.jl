

function getHarmonicFreqs(EoM!, bdys; kwargs...)
  _hbar = 1.0545718001391127e-34
  _e    = 1.6021766208e-19
  _amu  = 1.66053904e-27

  # Frequency conversion factor
  # first part is ASE style conversion
  # second part ( * 8065.61) is eV to cm-1
  c     = _hbar * 1e10 / sqrt(_e * _amu) * 8065.610420 

  #This code block is borrowed from opt 
  #Refer to optimizations.jl 
  x0         = prepX0(bdys)
  pars, mols = getPairs(bdys)
  vars       = optVars(mols, pars)
  im         = [i.m ^ -0.5 for i in bdys for j in 1:3]
  m          = im * im' #inverse mass scaling matrix

  function f(x, vars)
    G = zero(x)
    EoM!(nothing, G, x, vars)
    return G
  end

  H = jacobian(central_fdm(6,1), x -> f(x, vars), x0)[1]

  mH    = m .* H
  freqs = eigvals(mH)
  modes = eigvecs(mH)
  freqs = c * sqrt.(Complex.(freqs))

  return freqs, modes
end


