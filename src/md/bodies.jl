using StaticArrays

abstract type Particle end

abstract type Molecule end

#Atoms in simulation
struct Atom <: Particle
  r::SVector{3,Float64}
  v::SVector{3,Float64}
  m::Float64
  q::Float64
end

#CO molecule
struct CO <: Molecule
  C::Particle
  O::Particle
end
