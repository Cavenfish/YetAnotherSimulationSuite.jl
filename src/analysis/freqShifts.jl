
struct MolFreq
  Evib::Float64
  freq::Float64
  time::Float64
  r::Float64
  θ::Float64
end

function _getWave(bl)
  pks  = findPeaks(bl)
  i    = pks[1]
  j    = pks[2]
  wave = bl[i:j] .- mean(bl[i:j]) |> (x -> repeat(x, 1000))

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

    pks = findPeaks(I[1:n]; min=5e1)

    length(pks) > 0 || continue

    push!(freqs, MolFreq(Evib, v[pks[1]], tj.t[i], 0.0, 0.0))
  end

  freqs
end

function getAllFreqs(EoM, tj, dt)

  ind = let
    v = tj.v[102]
    v = [abs.(i) for i in v]

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

      pks = findPeaks(I[1:n]; min=5e1)

      length(pks) > 0 || continue

      sub = CoM(bdys[mol])
      d   = norm(main - sub)
      θ   = getAngleCO(bdys[ind], bdys[mol])
      push!(allFreqs, MolFreq(Evib, v[pks[1]], tj.t[i], d, θ))
    end
  end

  allFreqs
end

function getAllFreqs(EoM, bdys; ind=[0,1])

  allFreqs = MolFreq[]
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

    pks = findPeaks(I[1:n]; min=5e1)

    length(pks) > 0 || continue

    sub = CoM(bdys[mol])
    d   = norm(main - sub)
    θ   = getAngleCO(bdys[ind], bdys[mol])
    push!(allFreqs, MolFreq(Evib, v[pks[1]], 0.0, d, θ))
  end

  allFreqs
end

function getFvE(EoM, bdys, Erange)

  v,m = getHarmonicFreqs(EoM, bdys)

  freqs = []
  for E in Erange

    for bdy in bdys
      bdy.v .*= 0.0
    end

    vibExcite!(bdys, m[:,6], E)

    nve  = runNVE(EoM, (0, 50fs), 0.01fs, bdys; save="sparse") |> processDynamics
    # r    = [j[mol] for j in nve.r]
    bl   = [norm(j[2]-j[1]) for j in nve.r]
    wave = _getWave(bl)

    n = div(length(wave), 2)
    I = abs.(fft(wave))
    v = fftfreq(length(I), 1e17) ./ 29979245800.0

    pks = findPeaks(I[1:n]; min=1e3)
    
    push!(freqs, v[pks[1]])
  end

  freqs
end

function getFreqCoupling(EoM, bdysIn, E)
  bdys = deepcopy(bdysIn)

  ind  = pickRandomMol(bdys, "bulk") |> (x -> findall(e -> e in x, bdys))
  v,m  = getHarmonicFreqs(EoM, bdys[ind])

  b4   = getAllFreqs(EoM, bdys, ind=ind)
  r    = [i.r for i in b4]
  ϕ    = [i.θ for i in b4] .* 180/pi
  v    = [i.freq for i in b4]

  vibExcite!(bdys[ind], m[:,6], E)

  b5   = getAllFreqs(EoM, bdys, ind=ind)
  v2   = [i.freq for i in b5]
  dv   = v2 .- v

  (r, ϕ, dv)
end

