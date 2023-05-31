
struct Simulation
  bdys::Vector
  time::Vector
  energy::Vector
end

function run_NVE(EoM, tspan, dt, algo, bdys; kwargs...)
  pos   = [i.r for i in bdys]
  vel   = [i.v for i in bdys]
  m     = [i.m for i in bdys]
  simul = Simulation(bdys, [], [])

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simul; kwargs...)
  solu  = solve(prob, algo, dt=dt)

  return solu
end


#TODO:
# -Make callback function for simulations
# -callback can save values into struct