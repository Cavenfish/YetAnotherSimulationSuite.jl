# using OrdinaryDiffEq
# using StaticArrays

struct Simulation
  bdys::Vector
  time::Float64
  energy::Float64
end

function run_NVE(EoM, tspan, dt, algo, bdys)
  pos   = [i.r for i in bdys]
  vel   = [i.v for i in bdys]
  m     = [i.m for i in bdys]

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, m)
  solu  = solve(prob, algo, dt=dt)

  return solu
end


#TODO:
# -Make callback function for simulations
# -callback can save values into struct