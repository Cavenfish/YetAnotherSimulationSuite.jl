
function avgDFs(arr)
  ret = arr[1]

  for df in arr[2:end]
    ret .+= df
  end

  ret ./= length(arr)
  
  return ret
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

  return fig
end

function mkIsoJld(fname, db)

  function f(x)
    dfs = [jldopen(file)["df"] for file in x[2]]

    ret = x[1] => avgDFs(dfs)
    
    return ret
  end

  isos = map(f, collect(db))

  jldsave(fname; isos)
end

function pltIsotopes(db)

  function f(x)
    dfs = [jldopen(file)["df"] for file in x[2]]

    ret = x[1] => avgDFs(dfs)
    
    return ret
  end

  isos = map(f, collect(db))
  
  set_theme!(myLightTheme)

  fig = Figure()
  ax  = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Energy (eV)")

  for iso in isos
    x = iso[2].time ./ 1000 
    y = Float64.(iso[2].molVib)
    l = iso[1]

    lines!(ax, x, y, label=l)
  end

  axislegend(ax)

  return fig
end

function pltIsotopes(jldFile::String)
  isos = jldopen(jldFile)["isos"]

  set_theme!(myLightTheme)

  fig = Figure()
  ax  = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Energy (eV)")

  for iso in isos
    x = iso[2].time[102:end] ./ 1000 
    y = Float64.(iso[2].molVib[102:end])
    l = iso[1]

    lines!(ax, x, y, label=l)
  end

  axislegend(ax, position=:lb)

  return fig
end

function pltGeneralDisp(toPlt)

  set_theme!(myLightTheme)

  fig = Figure()
  ax  = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Energy (eV)")

  for p in toPlt
    
    l  = p[1]
    df = p[2]

    lines!(ax, df.time ./ 1000, df.molVib, label=l)
  end

  fig
end

function pltDoubleDisp(toPlt)

  set_theme!(myLightTheme)

  fig = Figure(size=(1,4))
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