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

struct NVE{AM<:AbstractMatrix}
  lattice::AM
end

"""
    NVE()

Construct a default NVE ensemble with zero lattice.
"""
NVE() = NVE(@SMatrix zeros(3,3))

"""
    NVE(lat::AbstractMatrix)

Construct an NVE ensemble from a lattice matrix.
"""
NVE(lat::AbstractMatrix) = SMatrix{size(lat)...}(lat) |> NVE

"""
    NVE(cell::MyCell)

Construct an NVE ensemble from a MyCell object.
"""
NVE(cell::MyCell) = NVE(cell.lattice)

struct NpT{B, AM<:AbstractMatrix, T<:MyThermostat}
  lattice::AM
  barostat::B
  thermostat::T
end

"""
    NVT{D,T,F}

Structure for NVT (canonical) ensemble.

# Fields
- `lattice`: Lattice matrix.
- `thermostat`: Thermostat object.
"""
struct NVT{AM<:AbstractMatrix, T<:MyThermostat}
  lattice::AM
  thermostat::T
end

"""
    NVT(thermostat::MyThermostat)

Construct an NVT ensemble with zero lattice and given thermostat.
"""
NVT(thermostat::MyThermostat) = NVT(zeros(3,3), thermostat)

"""
    NVT(cell::MyCell, thermostat::MyThermostat)

Construct an NVT ensemble from a MyCell and thermostat.
"""
NVT(cell::MyCell, thermostat::MyThermostat) = NVT(cell.lattice, thermostat)

"""
    NVT(lat::AbstractMatrix, thermostat::MyThermostat)

Construct an NVT ensemble from a lattice and thermostat.
"""
function NVT(lat::AbstractMatrix, thermostat::MyThermostat)
  l = SMatrix{size(lat)...}(lat)
  NVT(l, thermostat)
end

"""
    Dynamics{T,D,B,P,PV,I,F,S}

Structure holding all MD simulation variables.

# Fields
- `m`: Masses.
- `s`: Symbols.
- `mols`: Molecule indices.
- `temp`: Temperatures.
- `energy`: Energies.
- `forces`: Forces.
- `potVars`: Potential variables.
- `PBC`: Periodic boundary conditions.
- `NC`: Neighbor counts.
- `ensemble`: Ensemble object.
"""
struct Dynamics{T,D,B, PV<:PotVars, I<:Int, F<:AbstractFloat, S<:AbstractString}
  m::Vector{F}
  s::Vector{S}
  mols::Vector{Vector{I}}
  temp::Vector{F}
  energy::Vector{F}
  forces::Vector{Vector{SVector{D, F}}}
  potVars::PV
  PBC::Vector{B}
  NC::Vector{I}
  ensemble::T
end

"""
    singleRun(calc, vel, pos, tspan, simu, algo, dt; kwargs...)

Run a single MD simulation segment.

# Arguments
- `calc`: Calculator object.
- `vel`: Initial velocities.
- `pos`: Initial positions.
- `tspan`: Time span tuple.
- `simu`: Dynamics object.
- `algo`: ODE solver algorithm.
- `dt`: Time step.
- `kwargs`: Additional keyword arguments.

# Returns
- ODE solution object.
"""
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

"""
    doRun(calc, vel, pos, tspan, simu, algo, dt, split; kwargs...)

Run an MD simulation, optionally splitting into segments.

# Arguments
- `calc`: Calculator object.
- `vel`: Initial velocities.
- `pos`: Initial positions.
- `tspan`: Time span tuple.
- `simu`: Dynamics object.
- `algo`: ODE solver algorithm.
- `dt`: Time step.
- `split`: Number of segments.
- `kwargs`: Additional keyword arguments.

# Returns
- Processed trajectory or solution.
"""
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

"""
    run(calc, bdys, tspan, dt, ensemble; algo=VelocityVerlet(), split=1, kwargs...)

Run an MD simulation for a set of atoms.

# Arguments
- `calc`: Calculator object.
- `bdys`: Vector of MyAtoms.
- `tspan`: Time span tuple.
- `dt`: Time step.
- `ensemble`: Ensemble object.
- `algo`: ODE solver algorithm (default: VelocityVerlet()).
- `split`: Number of segments (default: 1).
- `kwargs`: Additional keyword arguments.

# Returns
- Processed trajectory or solution.
"""
function Base.run(
  calc::MyCalc, bdys::Vector{MyAtoms}, tspan::Tuple{Float64, Float64},
  dt::Float64, ensemble::T; algo=VelocityVerlet(), split=1, kwargs...
) where T <: Union{NVE,NpT,NVT}

  NC      = [0,0,0]
  PBC     = repeat([false], 3)
  symbols = [i.s for i in bdys]
  mas     = [i.m for i in bdys]
  potVars = calc.b(bdys)
  mols    = getMols(bdys, 1.2)
  pos     = [SVector{3}(i.r) for i in bdys]
  vel     = [SVector{3}(i.v) for i in bdys]

  simu = Dynamics(
    mas, symbols, mols, Float64[], Float64[], 
    Vector{SVector{3, Float64}}[], potVars, PBC, NC, ensemble
  )

  doRun(calc, vel, pos, tspan, simu, algo, dt, split; kwargs...)
end

"""
    run(calc, cell, tspan, dt, ensemble; algo=VelocityVerlet(), split=1, kwargs...)

Run an MD simulation for a cell.

# Arguments
- `calc`: Calculator object.
- `cell`: MyCell object.
- `tspan`: Time span tuple.
- `dt`: Time step.
- `ensemble`: Ensemble object.
- `algo`: ODE solver algorithm (default: VelocityVerlet()).
- `split`: Number of segments (default: 1).
- `kwargs`: Additional keyword arguments.

# Returns
- Processed trajectory or solution.
"""
function Base.run(
  calc::MyCalc, cell::MyCell, tspan::Tuple{Float64, Float64},
  dt::Float64, ensemble::T; algo=VelocityVerlet(), split=1, kwargs...
) where T <: Union{NVE,NpT,NVT}

  potVars = calc.b(cell)
  mols    = getMols(cell, 1.2)
  pos     = [SVector{3}(i) for i in getPos(cell)]
  vel     = [SVector{3}(i) for i in cell.velocity]

  simu = Dynamics(
    cell.masses, cell.symbols, mols, Float64[], Float64[], 
    Vector{SVector{3, Float64}}[], potVars, cell.PBC, cell.NC, ensemble
  )

  doRun(calc, vel, pos, tspan, simu, algo, dt, split; kwargs...)
end