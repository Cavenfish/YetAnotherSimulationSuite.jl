
#Two channel exponential decay
twoChannel(t, p) = p[1] * (exp.(- t / p[2]) + exp.(-t / p[3]))

#Exponential decay
expDecay(t, p) = p[1] * exp.(-t / p[2]) .+ p[3]

#Fit self dissipative cases
function selfDisp(df)

  t = df.time[102:end] ./ 1000 .- 10
  y = df.molVib[102:end]
  i = findall(e -> e > 0.055, y)
  j = findall(e -> e < 0.055, y)

  ifit = curve_fit(expDecay, t[i], y[i], [0.4, 500, 0.0])
  jfit = curve_fit(expDecay, t[j], y[j], [0.1, 50, 0.01])

  (ifit.param[2], jfit.param[2])
end

  