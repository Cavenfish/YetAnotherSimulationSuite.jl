
struct Thermostat{V,A, F<:Float64} <: MyThermostat
  T::F
  vars::V
  act!::A
end

function Thermostat(temp::Quantity, act!::Function; vars=nothing)
  T = uconvert(u"K", temp) |> ustrip
  
  Thermostat(T, vars, act!)
end

function getTemp(m, v, kB)
  # Leave out the 1/2 to get 2Ekin for T calc
  N    = length(m)
  Nf   = 3N - 3
  v2   = [i'i for i in v]
  Ekin = sum(m .* v2)
  
  Ekin / (kB * Nf)
end
