
function getDesorbEvents(files; dist=8.0)

  for file in files
    tj  = jldopen(file)["traj"]

    N   = length(tj.t)

    des = [] 

    for i in 1:N
      pos = tj.r[i]

      d   = [[norm(k - j) for j in pos] for k in pos]
      
      tmp = findall.(e -> e < 2, d)
      deleteat!.(d, tmp)

      num = findall(e -> minimum(e) > dist, d)

      push!(des, num)
    end

    fin = unique(vcat(des...))

    println("$file : ")
    println()
    if length(fin) > 0
      printstyled("  Desorbed Molecules: $(length(fin))"; bold=true, color=:green)
    else
      println("  Desorbed Molecules: $(length(fin))")
    end
  end
end