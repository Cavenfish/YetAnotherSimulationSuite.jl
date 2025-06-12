
function translate!(bdys::Vector{MyAtoms}, v::Vector{Float64})
  for i in bdys
    i.r .+= v
  end
end

function rotate!(bdys::Vector{MyAtoms}, R::Matrix)
  for i in bdys
    i.r .= R * i.r
  end
end

function rotate!(bdys::Vector{MyAtoms}, R::Matrix, about::Vector{Float64})
  translate!(bdys, -about)
  rotate!(bdys, R)
  translate!(bdys, about)
end

function rotate!(bdys::Vector{MyAtoms}, abc::Tuple{3,Float64})
  R = getRotationMatrix(abc...)
  rotate!(bdys, R)
end

function rotate!(bdys::Vector{MyAtoms}, abc::Tuple{3,F}, about::Vector{F}) where F<:AbstractFloat
  R = getRotationMatrix(abc...)
  rotate!(bdys, R, about)
end

function getRotationMatrix(α::F, β::F, γ::F) where F <: AbstractFloat
  Rz = [cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
  Ry = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)]
  Rx = [1 0 0; 0 cos(γ) -sin(γ); 0 sin(γ) cos(γ)]

  Rz*Ry*Rx
end

function randRotate!(bdys::Vector{MyAtoms})
  α = rand(-pi:1e-10:pi)
  γ = rand(-pi:1e-10:pi)
  β = rand(0:1e-10:pi)
  R = getRotationMatrix(α, β, γ)

  rotate!(bdys, R)
end

function randRotate!(bdys::Vector{MyAtoms}, about::Vector{Float64})
  α = rand(-pi:1e-10:pi)
  γ = rand(-pi:1e-10:pi)
  β = rand(0:1e-10:pi)
  R = getRotationMatrix(α, β, γ)

  rotate!(bdys, R, about)
end