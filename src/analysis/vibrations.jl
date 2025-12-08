"""
    getHarmonicFreqs(calc::MyCalc, obj::Union{MyCell, Vector{MyAtoms}})

Compute harmonic vibrational frequencies and modes for a cell or molecule.

# Arguments
- `calc`: Calculator object (`MyCalc`).
- `obj`: `MyCell` or vector of `MyAtoms`.

# Returns
- Tuple: (vector of frequencies, matrix of modes).
"""
function getHarmonicFreqs(calc::MyCalc, obj::Union{MyCell, Vector{MyAtoms}})
  _hbar = 1.0545718001391127e-34
  _e    = 1.6021766208e-19
  _amu  = 1.66053904e-27

  # Frequency conversion factor
  # first part is ASE style conversion
  # second part ( * 8065.61) is eV to cm-1
  c     = _hbar * 1e10 / sqrt(_e * _amu) * 8065.610420 

  # Prepare vars
  x0, vars = prep4pot(calc.b, obj)
  n        = length(x0)
  H        = zeros(n,n)
  cache    = JacobianCache(x0)

  # Get mass scaling matrix
  m = if isa(obj, MyCell)
    im = [i^-0.5 for i in obj.masses for j = 1:3]
    im * im'
  else
    im = [i.m^-0.5 for i in obj for j = 1:3]
    im * im'
  end

  # Barrier function
  function f(dx::Vector{Float64}, x::Vector{Float64})
    fg!(nothing, dx, x, vars, calc)
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

"""
    animateMode(bdys::Vector{MyAtoms}, mode, fileName; c=1.0)

Animate a vibrational mode and write the trajectory to a file.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
- `mode`: Mode vector to animate.
- `fileName`: Output file name.
- `c`: (Optional) Amplitude scaling factor (default: 1.0).
"""
function animateMode(bdys::Vector{MyAtoms}, mode, fileName; c=1.0)
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
        x,y,z = bdys[k].r .+ (sin(j) .* m[k] .* c)

        println(f, "$s   $x   $y   $z")
      end
    end
  end
  close(f)
end

"""
    getModePES(EoM, bdys, mode; range=collect(-1:0.001:2))

Compute the potential energy surface (PES) along a vibrational mode.

# Arguments
- `calc`: Calculator object (`MyCalc`).
- `obj`: `MyCell` or vector of `MyAtoms`.
- `mode`: Mode vector.
- `range`: (Optional) Range of displacements (default: -1:0.001:2).

# Returns
- Tuple: (`range`, energy values).
"""
function getModePES(calc::MyCalc, obj::Union{MyCell, Vector{MyAtoms}}, mode; range=collect(-1:0.001:2))

  x0, vars = prep4pot(calc.b, obj)

  y = zero(range)

  for (i,c) in enumerate(range)
    
    j = @. x0 + c * mode 

    y[i] = fg!(true, nothing, j, vars, calc)

  end

  (range, y)
end