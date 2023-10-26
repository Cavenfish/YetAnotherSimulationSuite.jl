
struct NVEsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
  m::Vector
end

struct NVTsimu
  bdys::Vector
  pars::Vector
  mols::Vector
  energy::Vector
  forces::Vector
  temp::Vector
  m::Vector
  thermostat!::Function
  thermoInps
end

function runNVE(EoM, tspan, dt, bdys; kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]   

  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [], mas)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)

  return solu
end

function runNVE(EoM, t1, t2, dt, bdys; kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]   
  tspan = (t1,t2)

  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [], mas)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false, save_start=false)

  return solu
end

function runNVT(EoM, tspan, dt, bdys, thermostat, thermoInps; kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]

  pars, mols = getPairs(bdys)
  simu       = NVTsimu(bdys, pars, mols, [], [], [], mas, thermostat, thermoInps)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt, dense=false, calck=false)

  return solu
end

#TODO:
# - Make Simulation struct more robust
# - Consider changing how bdys struct works