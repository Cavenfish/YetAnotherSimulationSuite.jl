
struct Thermostat{V,A, F<:AbstractFloat} <: MyThermostat
  T::F
  vars::V
  act!::A
end

function Thermostat(T::AbstractFloat, act!::Function; vars=nothing)
  Thermostat(T, vars, act!)
end

function getTemp(m, v)
  # Leave out the 1/2 to get 2Ekin for T calc
  N    = length(m)
  Nf   = 3N - 3
  v2   = [i'i for i in v]
  Ekin = sum(m .* v2)
  
  Ekin / (kB * Nf)
end
