struct CVR <:ThermoVars
  tau::Float64
end

function CVR(T::F, tau::F) where F<:AbstractFloat
  Thermostat(T, CVR!, vars=CVR(tau))
end

function CVR!(a, v, m, Tsim, thermostat)
  N    = length(m)
  tau  = thermostat.vars.tau

  if Tsim == 0.0
    return
  end

  Nf = (3N) - 3
  n  = Normal(0.0,1.0)
  R1 = rand(n)
  if (Nf-1)%2 == 0
    alpha = (Nf - 1)/2
    G     = Gamma(alpha, 1)
    R2    = 2 .* rand(G)
  else
    alpha = (Nf - 2)/2
    G     = Gamma(alpha, 1)
    R2    = 2 .* rand(G) .+ (rand(n) .^2)
  end

  if tau > 0.1
    c1 = exp(-1/tau)
  else
    c1 = 0.0
  end

  v2    = [i'i for i in v]
  K     = sum(0.5 .* m .* v2)
  sigma = 0.5 * Nf * kB * thermostat.T

  Knew   = @. ( 
                K + (1-c1) * (sigma * (R2 + R1^2) / Nf - K) + 
                2 * R1 * sqrt(K * sigma / Nf * (1-c1) * c1)
              )
  alpha2 = @. sqrt(Knew / K)

  @. v *= alpha2
end