
function spaghetti(files)

  N   = length(files)
  c   = 255
  k   = 220 / N
  f   = Figure()
  ax  = myAxis(f[1,1])  
  dfs = [jldopen(file)["df"] for file in files]

  for df in dfs

    x    = df.time ./ 1000
    y    = Float64.(df.molVib)

    lines!(ax, x, y, color=RGBf(c/255, 0, c/255))
    c -= k
  end

  x = dfs[1].time ./ 1000
  y = sum([df.molVib for df in dfs]) ./ length(dfs)

  lines!(ax, x, y, color=:gold)

  return f
end
    
