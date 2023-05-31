
struct Simulation
  bdys::Vector
  pars::Vector
  mols::Vector
  time::Vector
  energy::Vector
  forces::Vector
end

function run_NVE(EoM, tspan, dt, algo, bdys; kwargs...)
  pos   = [i.r for i in bdys]
  vel   = [i.v for i in bdys]
  m     = [i.m for i in bdys]

  pars, mols = getPairs(bdys)
  simul      = Simulation(bdys, pars, mols, [], [], [])

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simul; kwargs...)
  solu  = solve(prob, algo, dt=dt)

  return solu
end


#TODO:
# - Make Simulation struct more robust
# - Consider changing how bdys struct works