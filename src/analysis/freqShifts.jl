
struct MolFreq
  Evib::Float64
  freq::Float64
  time::Float64
  r::Float64
end

function getMolFreq(EoM, tj, dt)

  ind = let
    v = tj.v[102]

    i = findfirst(e -> e == maximum(v), v)

    isodd(i) ? [i, i+1] : [i-1, i]
  end 

  freqs = MolFreq[]
  for i in 1:dt:length(tj.t)

    bdys = getFrame(tj, i)
    Evib = getCOVibEnergy(byds[ind]; pot=EoM)

    nve  = runNVE(EoM, (0, 50fs), 0.1fs, bdys) |> (x -> processDynamics(x; dt=0.1fs))
    r    = [j[ind] for j in nve.r]
    bl   = [norm(j[2]-j[1]) for j in r]

    tmp  = findall(e -> isapprox(e, bl[1]; atol=3e-4), bl[2:end])

    wave = bl[1:tmp[2]] .- mean(bl[1:tmp[2]]) |> (x -> repeat(x, 500))

    n = div(length(wave), 2)
    I = abs.(fft(wave))
    v = fftfreq(length(I), 1e16) ./ 29979245800.0

    pks = findPeaks(I[1:n]; min=1e3)

    push!(freqs, MolFreq(Evib, I[pks[1]], tj.t[i], 0.0))
  end

  freqs
end