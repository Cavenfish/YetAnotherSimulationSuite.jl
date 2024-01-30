export getHarmonicFreqs, animateMode, getModePES

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

  fdm = central_fdm(7,1)
  H   = jacobian(fdm, x -> f(x, vars), x0)[1]

  mH    = m .* H
  mH    = Symmetric(mH) #Ensures eigvecs are not imaginary
  freqs = eigvals(mH)
  modes = eigvecs(mH)
  freqs = c * sqrt.(Complex.(freqs))

  return freqs, modes
end

function animateMode(bdys, mode, fileName)
  f = open(fileName, "w")
  N = length(bdys)

  m = Vector[]
  for i in 1:3:length(mode)
    push!(m, mode[i:i+2])
  end

  for i in 1:3
    for j in 0:pi/100:2*pi
      
      println(f, N)
      println(f, "Made by JMD")

      for k in 1:N
        
        s     = bdys[k].s
        x,y,z = bdys[k].r .+ sin(j) .* m[k]

        println(f, "$s   $x   $y   $z")
      end
    end
  end
  close(f)
end

function getModePES(EoM!, bdys, mode)

  x0, vars = prep4pot(bdys)

  x, y = [], []

  for i in -1:0.001:2
    
    j = @. x0 + i * mode 

    push!(x, i)
    push!(y, EoM!(true, nothing, j, vars))

  end
  return x, y
end

function getModeInteractionPES(EoM!, bdys, mode)

  x0, vars = prep4pot(bdys)
  molVar   = optVars(vars.mols, [])

  x, y = [], []

  for i in -1:0.001:2
    
    j = @. x0 + i * mode 

    E = EoM!(true, nothing, j, vars) - EoM!(true, nothing, j, molVar)

    push!(x, i)
    push!(y, E)

  end
  return x, y
end
    
