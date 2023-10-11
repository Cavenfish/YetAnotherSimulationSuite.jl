
struct DataBaseInfo
  experiment::String
  categories::String
end

struct exptInfo
  cluster::String
  molecule::String
  energy::Float64
end

function walk(db, group)
  ["$(group)/$(i)" for i in keys(db[group])]
end

function pull(db, group)
  [db["$(group)/$(i)"] for i in keys(db[group])]
end

function getMissing(db, group)
  i = keys(db[group]) |> (x -> parse.(Int64, x))
  j = collect(1:100)
  findall(e -> !(e in i), j)
end

function addAvgToDB!(db, group)
  avg = pull(db, group) |> avgDFs
  tmp = group[1:end-3]
  
  occursin("DF", group) ? key = tmp * "avgDF" : key = tmp * "avgVD"

  db[key] = avg
end

function makeAvgDB(file, db, p; kwargs...)
  
  jldopen(file, "a+"; kwargs...) do f

    avg = pull(db, p[1]) |> avgDFs

    f[p[2]] = avg
  end
end

function makeDataBase(file, p; kwargs...)

  jldopen(file, "a+"; kwargs...) do f 
    
    group = try 
              f[p[1]] 
            catch e 
              JLD2.Group(f, p[1]) 
            end 
    
    for item in p[2]
      k = split(item, "_")[end] |> (x -> replace(x, ".jld2" => ""))
      
      if k in keys(group)
        continue
      end

      if occursin("DF", p[1])
        df       = jldopen(item)["df"]
        group[k] = df
      elseif occursin("VD", p[1])
        vd       = jldopen(item)["vd"]
        group[k] = vd
      elseif occursin("RE", p[1])
        re       = jldopen(item)["re"]
        group[k] = re
      else
        error("Not VD or DF, so idk what to do")
      end
    end
  end
end

#"isotopes/13co-amCO/DFs" => glob("./gfs/isotopes/13co-amCO*DF*")

function makeDataBase(dbName::String, files::Vector{String}; kwargs...)

  jldopen(dbName, "a+"; kwargs...) do f 

    for file in files

      if "TJ" in file
        continue
      end

      headGroup = split(file, '/')[end-1]
      tailGroup = "DF" in file ? "DFs" : "VDs"
      k         = split(file, "_")[end] |> (x -> replace(x, ".jld2" => ""))

      if headGroup == "isotopes"

        midGroup = split(file, '/')[end] |> (x -> split(x, '_')[1])

        if midGroup == "surf"
        end #surf

        key = "$(headGroup)/$(midGroup)/$(tailGroup)/$(k)"
        tmp = "DF" in file ? "df" : "vd"

        f[key] = jldopen(file)[tmp]

      end #isotopes

    end # for file in files

  end #close JLD2 file
end