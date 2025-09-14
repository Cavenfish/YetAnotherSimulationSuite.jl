struct CVR <:ThermoVars
  tau::Float64
  kB::Float64
end

function CVR(T::Quantity, tau::Quantity, calc::MyCalc)
  t = uconvert(calc.time_unit, tau) |> ustrip

  Thermostat(T, CVR!, vars=CVR(t, calc.kB))
end

function CVR!(a, v, m, Tsim, thermostat)
  N   = length(m)
  tau = thermostat.vars.tau
  kB  = thermostat.vars.kB

  Tsim == 0.0 && return

  Nf = (3N) - 3
  n  = Normal(0.0,1.0)
  R1 = rand(n)
  
  R2 = if (Nf-1)%2 == 0
    alpha = (Nf - 1)/2
    G     = Gamma(alpha, 1)
    
    2 .* rand(G)
  else
    alpha = (Nf - 2)/2
    G     = Gamma(alpha, 1)
    
    2 .* rand(G) .+ (rand(n) .^2)
  end

  tau > 0.1 ? c1 = exp(-1/tau) : c1 = 0.0

  K     = 0.5 * Nf * kB * Tsim
  sigma = 0.5 * Nf * kB * thermostat.T

  Knew   = @. ( 
                K + (1-c1) * (sigma * (R2 + R1^2) / Nf - K) + 
                2 * R1 * sqrt(K * sigma / Nf * (1-c1) * c1)
              )
  alpha2 = @. sqrt(Knew / K)

  @. v *= alpha2
end