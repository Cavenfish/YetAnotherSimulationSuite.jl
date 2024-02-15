
function avgDFs(arr)
  ret = deepcopy(arr[1])

  for df in arr[2:end]
    ret .+= df
  end

  ret ./= length(arr)
  
  ret
end

function spaghetti(files)

  set_theme!(myLightTheme)

  N   = length(files)
  c   = 255
  k   = 220 / N
  fig = Figure()
  ax  = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Energy (eV)")
  dfs = [jldopen(file)["df"] for file in files]

  for df in dfs

    x    = df.time[102:end] ./ 1000
    y    = Float64.(df.molVib[102:end])

    lines!(ax, x, y, color=RGBf(c/255, 0, c/255))
    c -= k
  end

  x = dfs[1].time[102:end] ./ 1000
  y = sum([df.molVib[102:end] for df in dfs]) ./ length(dfs)

  lines!(ax, x, y, color=:gold, label="Average")

  axislegend(ax)

  fig
end

function pltGeneralDisp(toPlt)

  set_theme!(myLightTheme)

  fig = Figure()
  ax  = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Vibrational Energy (eV)")

  for p in toPlt
    
    l  = p[1]
    df = p[2]

    lines!(ax, df.time[102:end] ./ 1000, df.molVib[102:end], label=l)
  end

  fig
end

function pltDispAndTemp(toPlt)

  set_theme!(myLightTheme)

  fig = Figure(resolution=(800,700))
  gl  = GridLayout(fig[1, 1])

  ax1 = Axis(gl[1:6, 1], xlabel="Time (ps)", ylabel="Energy (eV)")
  ax2 = Axis(gl[7:8, 1], xlabel="Time (ps)", ylabel="Temp. (K)")

  hidexdecorations!(ax1, grid=false)

  for p in toPlt
    
    l  = p[1]
    df = p[2]

    lines!(ax1, df.time ./ 1e3, df.molVib, label=l)
    lines!(ax2, df.time ./ 1e3, df.temp)
  end

  axislegend(ax1)
  rowgap!(gl, 5)

  (fig, gl)
end

function pltDoubleDisp(toPlt)

  set_theme!(myLightTheme)

  fig = Figure(resolution=(1600,600))
  ax1 = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Energy (eV)")
  ax2 = Axis(fig[1,2], xlabel="Time (ps)", ylabel="Energy (eV)")

  for p in toPlt[1]
    
    l  = p[1]
    df = p[2]

    lines!(ax1, df.time ./ 1000, df.molVib, label=l)
  end

  for p in toPlt[2]
    
    l  = p[1]
    df = p[2]

    lines!(ax2, df.time ./ 1000, df.molVib, label=l)
  end

  fig
end

function pltGeneralVDOS(toPlt)

  set_theme!(myLightTheme)

  fig = Figure()
  ax  = Axis(fig[1,1], xlabel=L"Frequency (cm$^{-1]}$)",
                       ylabel="Intensity (arb.)")

  for p in toPlt
    
    l  = p[1]
    df = p[2]

    lines!(ax, df.v, df."1", label=l)
  end

  fig
end

function pltStackedVDOS(toPlt)

  set_theme!(myLightTheme)

  fig = Figure(resolution=(800, 200*length(toPlt)))
  gl  = GridLayout(fig[1, 1])

  for i in 1:length(toPlt)
    p   = toPlt[i]
    ax1 = Axis(gl[i, 1], xlabel=L"Frequency (cm$^{-1}$)", ylabel="VDOS (arb.)")
    ax2 = Axis(gl[i, 2], xlabel=L"Frequency (cm$^{-1}$)", ylabel="VDOS (arb.)")

    l  = p[1]
    df = p[2]
    y  = df."1" ./ maximum(df."1") .* 100

    lines!(ax1, df.v, y)
    xlims!(ax1, -5, 150)
    ylims!(ax1, -5, 100)

    lines!(ax2, df.v, y)
    xlims!(ax2, 2100, 2250)
    ylims!(ax2, -5, 100)


    # text!( ax, 1000, 50, text=l)

    hideydecorations!(ax2, grid=false)

    i < length(toPlt) && hidexdecorations!(ax1, grid=false)
    i < length(toPlt) && hidexdecorations!(ax2, grid=false)
  end

  rowgap!(gl, 0)

  (fig, gl)
end

function pltStackedRDF(toPlt)

  set_theme!(myLightTheme)

  fig = Figure(resolution=(600, 200*length(toPlt)))
  gl  = GridLayout(fig[1, 1])

  for i in 1:length(toPlt)
    p  = toPlt[i]
    ax = Axis(gl[i, 1], xlabel=L"Distance ($\AA$)", ylabel="g(r) (arb.)")

    l    = p[1]
    r, g = p[2]
    y    = g ./ maximum(g)

    lines!(ax, r, y)
    text!( ax, 20, 0.5, text=l, bold=true)
    xlims!(ax, 0, 30)
    ylims!(ax, -0.1, 1.1)

    i < length(toPlt) && hidexdecorations!(ax, grid=false)
  end

  rowgap!(gl, 0)

  (fig, gl)
end

function pltAvgEng(toPlt)

  set_theme!(myLightTheme)

  fig = Figure()
  ax  = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Translation Energy per Molecule (K)")

  for i in 1:length(toPlt)
    p = toPlt[i]

    l  = p[1]
    df = p[2]
    # tm = savGol(df.avgTra .- df.avgTra[1], 101, 5)
    x  = df.time[102:end] ./ 1000
    y  = (df.avgTra[102:end] .- df.avgTra[1]) .* 11600 .* 2 
    lines!(ax, x, y, label=l)
  end

  fig
end

function pltFreqShift(k)

  set_theme!(myLightTheme)

  fig = Figure()
  ax  = Axis(fig[1,1], xlabel="Energy (eV)", ylabel=L"Frequency (cm$^{-1}$)")

  contourf!(k.x, k.y, k.density, colormap=:acton)

  fig
end

struct TauVsDeltaNu
  label::String
  τ::Vector{Float64}
  τe::Vector{Float64}
  dv::Vector{Float64}
end

function getTauVsDeltaNu(df, vd, label; N=300)
 
  E     = df.molVib[102:end]
  τ, τe = trackTau(df, N)

  dv  = let
    i   = findall(e ->  1500< e < 3000, vd.v)
    v   = findPeaks(vd."1"[i], min=50, max=500)  |> (x -> vd.v[i[x[end]]])
    vs  = findPeaks(vd."1"[i], min=500) |> (x -> vd.v[i[x[1]]])
    v0  = vs / sqrt((11.23-E[1]) / 11.23)
    vp  = freqShiftMorse.([v0], [11.23], E)
    
    vp .- v
  end

  TauVsDeltaNu(label, τ, τe, dv[1:length(τ)])
end

function pltTauVsDeltaNu(toPlt)

  set_theme!(myLightTheme)

  fig = Figure(size=(800,600))
  ax  = Axis(fig[1,1], xlabel=L"$\Delta \nu$ (cm$^{-1}$)", ylabel=L"$\tau$ (ps)")

  for obj in toPlt
    l = convert(typeof(L""), obj.label)
    scatter!(ax, obj.dv, obj.τ, label=l)
  end

  fig
end

function pltTauVsDeltaNu(left, right)

  set_theme!(myLightTheme)

  fig = Figure(size=(1650,700))

  ax1 = Axis(fig[1,1], 
    xlabel=L"$\Delta \nu$ (cm$^{-1}$)", 
    ylabel=L"$\tau$ (ps)",
    xlabelsize=30, ylabelsize=30)

  ax2 = Axis(fig[1,2], 
    xlabel=L"$\Delta \nu$ (cm$^{-1}$)", 
    xlabelsize=30)

  for obj in left
    l = convert(typeof(L""), obj.label)
    scatter!(ax1, obj.dv, obj.τ, label=l)
  end

  for obj in right
    l = convert(typeof(L""), obj.label)
    scatter!(ax2, obj.dv, obj.τ, label=l)
  end

  fig[2, 1:2] = Legend(fig, ax1, orientation=:horizontal, 
                       labelsize=28, markersize=1000)

  hideydecorations!(ax2, grid=false)

  xlims!(ax1, -110, 110)
  ylims!(ax1, 0, 3000)

  xlims!(ax2, -110, 110)
  ylims!(ax2, 0, 3000)
  

  fig
end

