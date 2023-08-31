
function makeDataBase(file, p)

  jldopen(file, "w") do f 
    group = JLD2.Group(p[1])
    
    for item in p[2]
      k = split(item, "_")[end] |> (x -> replace(x, ".jld2" => ""))
      
      if occursin("DF", p[1])
        df       = jldopen(item)["df"]
        group[k] = df
      elseif occursin("VD", p[1])
        vd       = jldopen(item)["vd"]
        group[k] = vd
      else
        error("Not VD or DF, so idk what to do")
      end
    end
  end
end
