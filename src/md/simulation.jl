
struct NVEsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
  save::String
  m::Vector
  potVars::PotVars
end

struct NVTsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
  temp::Vector
  save::String
  m::Vector
  potVars::PotVars
  thermostat!::Function
  thermoInps
end

function runNVE(EoM, tspan, dt, bdys; save="full", kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]
  
  potVars    = EoM(bdys)
  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [], save, mas, potVars)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)

  return solu
end

function runNVE(EoM, t1, t2, dt, bdys; save="full", kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]
  μ     = [zeros(3) for i in bdys]
  tspan = (t1,t2)

  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [], save, mas, μ)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false, save_start=false)

  return solu
end

function runNVT(EoM, tspan, dt, bdys, thermostat, thermoInps; save="full", kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]
  μ     = [zeros(3) for i in bdys]

  pars, mols = getPairs(bdys)
  simu       = NVTsimu(bdys, pars, mols, [], [], [], save, mas, μ, thermostat, thermoInps)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)

  return solu
end

#TODO:
# - Make Simulation struct more robust
# - Consider changing how bdys struct works