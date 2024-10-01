
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

struct TauVsDeltaNu
  label::String
  τ::Vector{Float64}
  τe::Vector{Float64}
  dv::Vector{Float64}
end

function trackTau(df, N)
  t,y = prep4fit(df)

  τ   = Float64[]
  τe  = Float64[]
  for i in 0:length(t)- N - 1
    fit = curve_fit(expDecay, t[1+i:N+i], y[1+i:N+i], [0.4, 500.0])
    err = standard_errors(fit)

    if err[2] / fit.param[2] > 0.05
      return τ, τe
    end

    push!(τ , fit.param[2])
    push!(τe, err[2])
  end

  τ, τe
end

function getTauVsDeltaNu(df, vd, label; N=300)
 
  E     = df.molVib[102:end]
  τ, τe = trackTau(df, N)

  dv  = let
    i   = findall(e ->  1500< e < 3000, vd.v)
    v   = findPeaks(vd."1"[i], min=10, max=500, width=10)  |> (x -> vd.v[i[x[end]]])
    vs  = findPeaks(vd."1"[i], min=500) |> (x -> vd.v[i[x[1]]])
    v0  = vs / sqrt((11.23-E[1]) / 11.23)
    vp  = freqShiftMorse.([v0], [11.23], E)
    
    vp .- v
  end

  TauVsDeltaNu(label, τ, τe, dv[1:length(τ)])
end