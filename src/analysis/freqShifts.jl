
struct MolFreq
  Evib::Float64
  freq::Float64
  time::Float64
  r::Float64
  θ::Float64
end

struct INM
  E::Float64
  ν::Float64
end

function getbl(nve, vec)
  foo(x) = [j for i in x for j in i]
  bar(y) = foo(y) |> (x -> dot(x, vec))

  [bar(nve.r[i]) for i = 1:length(nve.r)]
end

function getWave(bl)
  pks  = findTurningPoints(bl)
  i    = pks[1]
  j    = pks[3]
  wave = bl[i:j] .- mean(bl[i:j]) |> (x -> repeat(x, 1000))

  wave
end

function getINM(EoM, bdys, mol, eignvec, energies; t=50fs, dt=0.01fs)
  inms = INM[]

  for E in energies

    tmp  = deepcopy(bdys)
    vibExcite!(tmp[mol], eignvec, E)


    nve  = runNVE(EoM, (0, t), dt, tmp; save="sparse") |> processDynamics
    frms = getFrame.([nve], collect(1:length(nve.t)))
    x1   = [j for i in frms[1][mol] for j in i.r]
    
    bl = Float64[]
    for frm in frms[2:end]
      xi = [j for i in frm[mol] for j in i.r]
      dot(xi - x1, eignvec) |> (x -> push!(bl, x))
    end
    
    wave = try
      getWave(bl)
    catch err
      @warn "Skipping Energy: $(E)"
      continue
    end

    n = div(length(wave), 2)
    I = abs.(fft(wave))
    p = 1fs/dt |> round |> (x -> 1e15 * x)
    v = fftfreq(length(I), p) ./ 29979245800.0

    pks = findPeaks(I[1:n]; min=5e1)

    length(pks) > 0 || continue

    push!(inms, INM(E, v[pks[1]]))
  end

  inms
end

function getINM(EoM, mol, eignvec, tj; dt=1, nveTime=50fs, nveDt=0.001fs)
  inms = INM[]

  for i in 1:dt:length(tj.t)
    
    bdys = getFrame(tj, i)
    
    #NEED TO FIX
    E = getVibEnergy(bdys[mol], eignvec, pot=EoM)

    nve  = runNVE(EoM, (0, nveTime), nveDt, bdys; save="sparse") |> processDynamics
    frms = getFrame.([nve], collect(1:length(nve.t)))
    x1   = [j for i in frms[1] for j in i.r]
    
    bl = Float64[]
    for frm in frms
      xi = [j for i in frm for j in i.r]
      dot(xi - x1, eignvec) |> (x -> push!(bl, x))
    end
    
    wave = try
      getWave(bl)
    catch err
      @warn "Skipping timestep: $(tj.t[i])"
      continue
    end

    n = div(length(wave), 2)
    I = abs.(fft(wave))
    p = 1fs/nveDt |> round |> (x -> 1e15 * x)
    v = fftfreq(length(I), p) ./ 29979245800.0

    pks = findPeaks(I[1:n]; min=5e1)

    length(pks) > 0 || continue

    push!(inms, INM(E, v[pks[1]]))
  end

  inms
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
      getWave(bl)
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
        getWave(bl)
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
      getWave(bl)
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

function getFvE(EoM, bdys, Erange, vec)

  freqs = []
  for E in Erange

    for bdy in bdys
      bdy.v .*= 0.0
    end

    try
      vibExcite!(bdys, vec, E)

      nve  = runNVE(EoM, (0, 100fs), 0.01fs, bdys; save="sparse") |> processDynamics
      bl   = getbl(nve, vec)
      wave = getWave(bl)

      n = div(length(wave), 2)
      I = abs.(fft(wave))
      v = fftfreq(length(I), 1e17) ./ 29979245800.0

      pks = findPeaks(I[1:n]; min=1e3)
      
      push!(freqs, v[pks[1]])
    catch
      @warn "Had to skip energy E=$(E)"
      push!(freqs, NaN)
    end

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

