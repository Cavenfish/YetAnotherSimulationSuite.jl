
struct MolFreq
  Evib::Float64
  freq::Float64
  time::Float64
  r::Float64
end

function _getWave(bl)
  tmpi = findall(e -> isapprox(e, maximum(bl); atol=5e-5), bl)
  tmpj = findall(e -> 1000 < e - tmpi[1] < 2000, tmpi)

  i    = tmpi[1]
  j    = tmpi[tmpj[div(length(tmpj), 2)]]
  wave = bl[i:j] .- mean(bl[i:j]) |> (x -> repeat(x, 500))

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

    nve  = runNVE(EoM, (0, 50fs), 0.01fs, bdys; save="sparse") |> processDynamics
    r    = [j[ind] for j in nve.r]
    bl   = [norm(j[2]-j[1]) for j in r]

    wave = try
      _getWave(bl)
    catch err
      @warn "Skipping timestep: $(tj.t[i])"
      continue
    end

    n = div(length(wave), 2)
    I = abs.(fft(wave))
    v = fftfreq(length(I), 1e17) ./ 29979245800.0

    pks = findPeaks(I[1:n]; min=1e1)

    length(pks) > 0 || continue

    push!(freqs, MolFreq(Evib, v[pks[1]], tj.t[i], 0.0))
  end

  freqs
end

function getAllFreqs(EoM, tj, dt)

  ind = let
    v = tj.v[102]

    i = findfirst(e -> e == maximum(v), v)

    isodd(i) ? [i, i+1] : [i-1, i]
  end 

  allFreqs = MolFreq[]
  for i in 1:dt:length(tj.t)

    bdys = getFrame(tj, i)
    nve  = runNVE(EoM, (0, 50fs), 0.01fs, bdys; save="sparse") |> processDynamics
    mols = getMols(bdys)
    main = CoM(bdys[ind])

    for mol in mols
      Evib = getCOVibEnergy(bdys[mol]; pot=EoM)
      r    = [j[mol] for j in nve.r]
      bl   = [norm(j[2]-j[1]) for j in r]

      wave = try
        _getWave(bl)
      catch err
        @warn "Skipping timestep: $(tj.t[i])"
        continue
      end

      n = div(length(wave), 2)
      I = abs.(fft(wave))
      v = fftfreq(length(I), 1e17) ./ 29979245800.0

      pks = findPeaks(I[1:n]; min=1e1)

      length(pks) > 0 || continue

      sub = CoM(bdys[mol])
      d   = norm(main - sub)
      push!(allFreqs, MolFreq(Evib, v[pks[1]], tj.t[i], d))
    end
  end

  allFreqs
end