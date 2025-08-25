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

struct NpT{D,B, T<:MyThermostat, F<:AbstractFloat}
  lattice::MMatrix{D, D, F}
  barostat::B
  thermostat::T
end

struct NVT{D, T<:MyThermostat, F<:AbstractFloat}
  lattice::SMatrix{D, D, F}
  thermostat::T
end

NVT(thermostat::MyThermostat) = NVT(zeros(3,3), thermostat)
NVT(cell::MyCell, thermostat::MyThermostat) = NVT(cell.lattice, thermostat)

function NVT(lat::AbstractMatrix, thermostat::MyThermostat)
  l = SMatrix{size(lat)...}(lat)
  NVT(l, thermostat)
end

struct Dynamics{T,D,B,P, PV<:PotVars, I<:Int, F<:AbstractFloat, S<:AbstractString}
  m::Vector{F}
  s::Vector{S}
  pars::Vector{P}
  mols::Vector{Vector{I}}
  temp::Vector{F}
  energy::Vector{F}
  forces::Vector{Vector{MVector{D, F}}}
  potVars::PV
  PBC::Vector{B}
  NC::Vector{I}
  ensemble::T
end

function singleRun(
  calc::MyCalc, vel::T, pos::T, tspan::Tuple{Float64, Float64},
  simu::Dynamics, algo::A, dt::AbstractFloat; kwargs...
) where {T,A}

  prob  = SecondOrderODEProblem(
    (dv, v, u, p, t) -> dyn!(dv, v, u, p, t, calc), 
    vel, pos, tspan, simu; kwargs...
  )
  
  solve(prob, algo; dt=dt, dense=false, calck=false)
end

function doRun(
  calc::MyCalc, vel::T, pos::T, tspan::Tuple{Float64, Float64},
  simu::Dynamics, algo::A, dt::AbstractFloat, split::Int; kwargs...
) where {T,A}

  if split > 1
    t = tspan[1]
    Δ = (tspan[2] - tspan[1]) / split
    
    for i = 1:split
      span = (t, t+Δ)
      solu = singleRun(calc, vel, pos, span, simu, algo, dt; kwargs...)

      open("./$(i).tmp", "w") do io
        serialize(io, solu)
      end

      # Is this needed?
      @free solu

      t   += Δ
    end

    files = ["./$(i).tmp" for i = 1:split]
    return processTmpFiles(files; dt=dt)
  else
    solu = singleRun(calc, vel, pos, tspan, simu, algo, dt; kwargs...)
    return processDynamics(solu)
  end

end

function Base.run(
  calc::MyCalc, bdys::Vector{MyAtoms}, tspan::Tuple{Float64, Float64},
  dt::Float64, ensemble::T; algo=VelocityVerlet(), split=1, kwargs...
) where T <: Union{NVE,NpT,NVT}

  NC         = [0,0,0]
  PBC        = repeat([false], 3)
  symbols    = [i.s for i in bdys]
  mas        = [i.m for i in bdys]
  potVars    = calc.b(bdys)
  pars, mols = getPairs(bdys)
  pos        = [MVector{3}(i.r) for i in bdys]
  vel        = [MVector{3}(i.v) for i in bdys]

  simu = Dynamics(
    mas, symbols, pars, mols, Float64[], Float64[], 
    Vector{MVector{3, Float64}}[], potVars, PBC, NC, ensemble
  )

  doRun(calc, vel, pos, tspan, simu, algo, dt, split; kwargs...)
end

function Base.run(
  calc::MyCalc, cell::MyCell, tspan::Tuple{Float64, Float64},
  dt::Float64, ensemble::T; algo=VelocityVerlet(), split=1, kwargs...
) where T <: Union{NVE,NpT,NVT}

  potVars    = calc.b(cell)
  pars, mols = getPairs(cell)
  pos        = [MVector{3}(i) for i in getPos(cell)]
  vel        = [MVector{3}(i) for i in cell.velocity]

  simu = Dynamics(
    cell.masses, cell.symbols, pars, mols, Float64[], Float64[], 
    Vector{MVector{3, Float64}}[], potVars, cell.PBC, cell.NC, ensemble
  )

  doRun(calc, vel, pos, tspan, simu, algo, dt, split; kwargs...)
end