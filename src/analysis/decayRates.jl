
#Two channel exponential decay
twoChannel(t, p) = (p[1] * exp.(- t / p[2])) .+ (p[3] * exp.(-t / p[4])) .+ abs(p[5])

#Exponential decay
expDecay(t, p) = p[1] * exp.(-t / p[2]) #.+ abs(p[3])

#Exponential decay with sin function for random gains
expSinDecay(t, p) = p[1] * exp.(-t / p[2]) .+ (p[3] * sin.(p[4] * t))

#My custom decay-rate with frequency gap law included
function myDecay(t,p)
  τ =  p[1] * exp.(t / p[2])
  0.405 * exp.(-t ./ τ) .+ 0.01
 end

function prep4fit(df)
  t = df.time[102:end] ./ 1000 .- 10
  y = df.molVib[102:end]

  (t,y)
end

function trackTau(df; N=100)
  t,y = prep4fit(df)

  X   = Iterators.partition(t, N) |> collect
  Y   = Iterators.partition(y, N) |> collect

  τ   = []
  for i in 1:length(X)
    fit = curve_fit(expDecay, X[i], Y[i], [0.4, 500.0, 0.0])
    push!(τ, fit.param[2])
  end

  τ
end

function trackTau(df, N)
  t,y = prep4fit(df)

  τ   = []
  for i in 0:length(t)- N - 1
    fit = curve_fit(expDecay, t[1+i:N+i], y[1+i:N+i], [0.4, 500.0])
    push!(τ, fit.param[2])
  end

  τ
end

#Fit self dissipative cases
function selfDisp(df; thresh=0.055)

  t = df.time[102:end] ./ 1000 .- 10
  y = df.molVib[102:end]
  i = findall(e -> e > thresh, y)
  j = findall(e -> e < thresh, y)

  ifit = curve_fit(expDecay, t[i], y[i], [0.4, 500, 0.0])
  jfit = curve_fit(expDecay, t[j], y[j], [0.1, 50, 0.01])

  (ifit, jfit)
end

function nonSelfDisp(df)

  t = df.time[102:end] ./ 1000 .- 10
  y = df.molVib[102:end]

  afit = curve_fit(expDecay, t, y, [0.4, 500, 1500, ])

  afit
end