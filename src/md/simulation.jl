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

struct NVE end

struct NVT{C, TV<:ThermoVars}
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
  lattice::SMatrix{D, D, F}
  ensemble::T
end

function Base.run(
  EoM::Function, bdys::Vector{MyAtoms}, tspan::Tuple{Float64, Float64},
  dt::Float64, ensemble::T; kwargs...
) where T <: Union{NVE,NVT}
             
  NC         = [0,0,0]
  PBC        = repeat([false], 3)
  lattice    = @SMatrix zeros(3,3)
  symbols    = [i.s for i in bdys]
  mas        = [i.m for i in bdys]
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  pos        = [SVector{3}(i.r) for i in bdys]
  vel        = [SVector{3}(i.v) for i in bdys]

  simu = Dynamics(
    mas, symbols, pars, mols, Float64[], Float64[], [],
    potVars, PBC, NC, lattice, ensemble
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)
end

function runMD(EoM, cell::MyCell, tspan::Tuple{Float64, Float64}, dt::Float64; 
               save="full", thermostat=nothing, thermoinps=nothing, kwargs...)

  bdys       = makeBdys(cell)
  potVars    = EoM(cell)
  pars, mols = getPairs(bdys)
  pos        = [SVector{3}(i.r) for i in bdys]
  vel        = [SVector{3}(i.v) for i in bdys]

  if thermoinps != nothing && thermostat != nothing
    NVT = true
  else
    NVT = false
  end

  simu = Simulation(
    bdys, pars, mols, [], [], [],
    save, cell.PBC, cell.NC, cell.lattice, cell.masses, potVars,
    NVT, thermostat, thermoinps
  )

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  
  solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)
end