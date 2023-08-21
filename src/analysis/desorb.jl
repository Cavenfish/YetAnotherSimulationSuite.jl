
function getDesorbEvents(files)

  for file in files
    tj  = jldopen(file)["traj"]

    N   = length(tj.t)

    des = [] 

    for i in 1:N
      pos = tj.r[i]

      d   = [[norm(k - j) for j in pos] for k in pos]
      
      tmp = findall.(e -> e < 2, d)
      deleteat!.(d, tmp)

      num = findall(e -> minimum(e) > 8, d)

      push!.(des, num)
    end

    des = unique(des)

    println("$file : ")
    println("  Desorbed Molecules: $(length(des))")
  end
end