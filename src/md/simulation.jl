
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
  energy::Vector{Float64}
  forces::Vector
  m::Vector
  temp::Vector{Float64}
  thermostat!::Function
  thermoInps
end

function runNVE(EoM, tspan, dt, bdys; kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]
  mas   = [i.m for i in bdys]

  #Consider swapping to this format
  #Might be much nicer to work with
  #-------------
  # N = length(pos)
  # p = zeros((3,N))
  # v = zeros((3,N))
  # for i in 1:N
  #   @views p[:,i] = pos[i][:]
  #   @views v[:,i] = vel[i][:]
  # end    

  pars, mols = getPairs(bdys)
  simu       = NVEsimu(bdys, pars, mols, [], [], mas)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt)

  return solu
end

function runNVT(EoM, tspan, dt, bdys, thermostat, thermoInps; kwargs...)
  pos   = [SVector{3}(i.r) for i in bdys]
  vel   = [SVector{3}(i.v) for i in bdys]

  pars, mols = getPairs(bdys)
  simu       = NVTsimu(bdys, pars, mols, [], [], [], thermostat, thermoInps)

  prob  = SecondOrderODEProblem(EoM, vel, pos, tspan, simu; kwargs...)
  solu  = solve(prob, VelocityVerlet(), dt=dt)

  return solu
end

#TODO:
# - Make Simulation struct more robust
# - Consider changing how bdys struct works