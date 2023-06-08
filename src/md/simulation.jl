
struct NVEsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
end

struct NVTsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector{Float64}
  forces::Vector
  temp::Vector{Float64}
  thermostat!::Function
  thermoInps
end

function runNVE(EoM, tspan, dt, bdys; kwargs...)
  pos   = [i.r for i in bdys]
  vel   = [i.v for i in bdys]

  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [])

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt)

  return solu
end

function runNVT(EoM, tspan, dt, bdys, thermostat, thermoInps; kwargs...)
  pos   = [i.r for i in bdys]
  vel   = [i.v for i in bdys]

  pars, mols = getPairs(bdys)
  simu       = NVTsimu(bdys, pars, mols, [], [], [], thermostat, thermoInps)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt)

  return solu
end

#TODO:
# - Make Simulation struct more robust
# - Consider changing how bdys struct works