
struct DataBaseInfo
  experiment::String
  categories::String
end

struct exptInfo
  cluster::String
  molecule::String
  energy::Float64
end

function makeDataBase(file, p; kwargs...)

  jldopen(file, "a+"; kwargs...) do f 
    group = JLD2.Group(f, p[1])
    
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
