"""
I want to refactor this to resemble the Traj style. One struct with
universal properties, and multiple other structs with more specific 
properties. This will reduce unnecessary memory.

Notes
----------
The simulations are where performance matters, hence any and all performance
tips and tricks should be done here. Making r and v SVectors prior to the
dynamics run is one example. 
"""

struct NVE{D, F<:AbstractFloat}
  lattice::SMatrix{D, D, F}
end

NVE() = NVE(@SMatrix zeros(3,3))
NVE(lat::AbstractMatrix) = SMatrix{size(lat)...}(lat) |> NVE
NVE(cell::MyCell) = NVE(cell.lattice)

struct NpT{C,D, TV<:ThermoVars, BV<:BaroVars, F<:AbstractFloat}
  lattice::MMatrix{D, D, F}
  baroVars::BV
  barostat!::C
  thermoVars::TV
  thermostat!::C
end

# C is callable
struct NVT{C,D, TV<:ThermoVars, F<:AbstractFloat}
  lattice::SMatrix{D, D, F}
  thermoVars::TV
  thermostat!::C
end

struct Dynamics{T,D,B,P, PV<:PotVars, I<:Int, F<:AbstractFloat, S<:AbstractString}
  m::Vector{F}
  s::Vector{S}
  pars::Vector{P}
  mols::Vector{Vector{I}}
  temp::Vector{F}
  energy::Vector{F}
  forces::Vector{Vector{SVector{D, F}}}
  potVars::PV
  PBC::Vector{B}
  NC::Vector{I}
  ensemble::T
end

function Base.run(
  EoM::Function, bdys::Vector{MyAtoms}, tspan::Tuple{Float64, Float64},
  dt::Float64, ensemble::T; algo=VelocityVerlet(), kwargs...
) where T <: Union{NVE,NpT,NVT}

  NC         = [0,0,0]
  PBC        = repeat([false], 3)
  symbols    = [i.s for i in bdys]
  mas        = [i.m for i in bdys]
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  pos        = [SVector{3}(i.r) for i in bdys]
  vel        = [SVector{3}(i.v) for i in bdys]

  simu = Dynamics(
    mas, symbols, pars, mols, Float64[], Float64[], 
    Vector{SVector{3, Float64}}[], potVars, PBC, NC, ensemble
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, algo, dt=dt, dense=false, calck=false)
end

function Base.run(
  EoM::Function, cell::MyCell, tspan::Tuple{Float64, Float64},
  dt::Float64, ensemble::T; algo=VelocityVerlet(), kwargs...
) where T <: Union{NVE,NpT,NVT}

  potVars    = EoM(cell)
  pars, mols = getPairs(cell)
  pos        = [SVector{3}(i) for i in getPos(cell)]
  vel        = [SVector{3}(i) for i in cell.velocity]

  simu = Dynamics(
    cell.masses, cell.symbols, pars, mols, Float64[], Float64[], 
    Vector{SVector{3, Float64}}[], potVars, cell.PBC, cell.NC, ensemble
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, algo, dt=dt, dense=false, calck=false)
end