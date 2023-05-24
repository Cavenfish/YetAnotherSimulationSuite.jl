using OrdinaryDiffEq

function run_NVE(system, EoM, tspan, dt, algo)

  prob = SecondOrderODEProblem(EoM, vel, pos, tspan)
  solu = solve(prob, algo, dt)

end
