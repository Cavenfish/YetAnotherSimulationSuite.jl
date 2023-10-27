
struct MolFreq
  Evib::Float64
  freq::Float64
  time::Float64
  r::Float64
end

function _getWave(bl)
  aTol = bl[1] / 1e6
  tmp  = findall(e -> isapprox(e, bl[1]; atol=aTol), bl[10:end])

  while length(tmp) < 2
    println(aTol)
    aTol *= 1.2
    tmp   = findall(e -> isapprox(e, bl[1]; atol=aTol), bl[10:end])
  end
  
  wave = bl[1:tmp[2]+9] .- mean(bl[1:tmp[2]+9]) |> (x -> repeat(x, 500))

  wave
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
    Evib = getCOVibEnergy(bdys[ind]; pot=EoM)

    nve  = runNVE(EoM, (0, 50fs), 0.1fs, bdys) |> (x -> processDynamics(x; dt=0.1fs))
    r    = [j[ind] for j in nve.r]
    bl   = [norm(j[2]-j[1]) for j in r]

    wave = _getWave(bl)
    return wave

    n = div(length(wave), 2)
    I = abs.(fft(wave))
    v = fftfreq(length(I), 1e16) ./ 29979245800.0

    pks = findPeaks(I[1:n]; min=1e1)

    length(pks) > 0 || continue
    println(v[pks[1]])

    push!(freqs, MolFreq(Evib, v[pks[1]], tj.t[i], 0.0))
  end

  freqs
end