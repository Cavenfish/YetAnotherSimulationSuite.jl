"""
    translate!(bdys::Vector{MyAtoms}, v::Vector{Float64})

Translate all atoms in `bdys` by vector `v` in-place.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
- `v`: Translation vector.
"""
function translate!(bdys::Vector{MyAtoms}, v::Vector{Float64})
  for i in bdys
    i.r .+= v
  end
end

"""
    LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, R::Matrix)

Rotate all atoms in `bdys` by rotation matrix `R` in-place.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
- `R`: 3x3 rotation matrix.
"""
function LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, R::Matrix)
  for i in bdys
    i.r .= R * i.r
  end
end

"""
    LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, R::Matrix, about::Vector{Float64})

Rotate all atoms in `bdys` by rotation matrix `R` about a point `about` in-place.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
- `R`: 3x3 rotation matrix.
- `about`: Point to rotate about.
"""
function LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, R::Matrix, about::Vector{Float64})
  translate!(bdys, -about)
  rotate!(bdys, R)
  translate!(bdys, about)
end

"""
    LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, abc::Tuple{F,F,F}) where F <: AbstractFloat

Rotate all atoms in `bdys` by Euler angles `(α, β, γ)` in-place.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
- `abc`: Tuple of Euler angles `(α, β, γ)` in radians.
"""
function LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, abc::Tuple{F,F,F}) where F <: AbstractFloat
  R = getRotationMatrix(abc...)
  rotate!(bdys, R)
end

"""
    LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, abc::Tuple{F,F,F}, about::Vector{F}) where F<:AbstractFloat

Rotate all atoms in `bdys` by Euler angles `(α, β, γ)` about a point `about` in-place.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
- `abc`: Tuple of Euler angles `(α, β, γ)` in radians.
- `about`: Point to rotate about.
"""
function LinearAlgebra.rotate!(bdys::Vector{MyAtoms}, abc::Tuple{F,F,F}, about::Vector{F}) where F<:AbstractFloat
  R = getRotationMatrix(abc...)
  rotate!(bdys, R, about)
end

"""
    getRotationMatrix(α::F, β::F, γ::F) where F <: AbstractFloat

Construct a rotation matrix from Euler angles `(α, β, γ)`.

# Arguments
- `α`, `β`, `γ`: Euler angles in radians.

# Returns
- 3x3 rotation matrix.
"""
function getRotationMatrix(α::F, β::F, γ::F) where F <: AbstractFloat
  Rz = [cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
  Ry = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)]
  Rx = [1 0 0; 0 cos(γ) -sin(γ); 0 sin(γ) cos(γ)]

  Rz*Ry*Rx
end

"""
    randRotate!(bdys::Vector{MyAtoms})

Apply a random rotation to all atoms in `bdys` in-place.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
"""
function randRotate!(bdys::Vector{MyAtoms})
  α = rand(-pi:1e-10:pi)
  γ = rand(-pi:1e-10:pi)
  β = rand(0:1e-10:pi)
  R = getRotationMatrix(α, β, γ)

  rotate!(bdys, R)
end

"""
    randRotate!(bdys::Vector{MyAtoms}, about::Vector{Float64})

Apply a random rotation about a point `about` to all atoms in `bdys` in-place.

# Arguments
- `bdys`: Vector of `MyAtoms` objects.
- `about`: Point to rotate about.
"""
function randRotate!(bdys::Vector{MyAtoms}, about::Vector{Float64})
  α = rand(-pi:1e-10:pi)
  γ = rand(-pi:1e-10:pi)
  β = rand(0:1e-10:pi)
  R = getRotationMatrix(α, β, γ)

  rotate!(bdys, R, about)
end