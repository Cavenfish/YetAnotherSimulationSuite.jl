

function getGrad(EoM!, bdys; kwargs...)
  #This code block is borrowed from opt 
  #Refer to optimizations.jl 
  x0         = prep_x0(bdys)
  pars, mols = getPairs(bdys)
  vars       = optVars(mols, pars)

  # Pre-allocate gradient
  G = zero(x0)
  EoM!(nothing, G, x0, vars)

  return G
end
