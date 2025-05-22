
function getHarmonicFreqs(EoM, bdys; kwargs...)
  _hbar = 1.0545718001391127e-34
  _e    = 1.6021766208e-19
  _amu  = 1.66053904e-27

  # Frequency conversion factor
  # first part is ASE style conversion
  # second part ( * 8065.61) is eV to cm-1
  c     = _hbar * 1e10 / sqrt(_e * _amu) * 8065.610420 

  # Prepare vars
  x0, vars = prep4pot(EoM, bdys)
  im       = [i.m ^ -0.5 for i in bdys for j in 1:3]
  m        = im * im' #inverse mass scaling matrix
  n        = length(x0)
  H        = zeros(n,n)
  cache    = JacobianCache(x0)

  # Barrier function
  function f(dx::Vector{Float64}, x::Vector{Float64})
    EoM(nothing, dx, x, vars)
  end

  # Calculate Hessian
  finite_difference_jacobian!(H, f, x0, cache)

  mH    = m .* H
  mH    = Symmetric(mH) #Ensures eigvecs are not imaginary
  freqs = eigvals(mH)
  modes = eigvecs(mH)
  freqs = c * sqrt.(Complex.(freqs))

  (freqs, modes)
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

function getModePES(EoM, bdys, mode; range=collect(-1:0.001:2))

  x0, vars = prep4pot(EoM, bdys)

  x, y = [], []

  for i in range
    
    j = @. x0 + i * mode 

    push!(x, i)
    push!(y, EoM(true, nothing, j, vars))

  end
  return x, y
end

function getModeInteractionPES(EoM!, bdys, mode; range=collect(-1:0.001:2))

  x0, vars = prep4pot(bdys)
  molVar   = optVars(vars.mols, [])

  x, y = [], []

  for i in range
    
    j = @. x0 + i * mode 

    E = EoM!(true, nothing, j, vars) - EoM!(true, nothing, j, molVar)

    push!(x, i)
    push!(y, E)

  end
  return x, y
end
    
