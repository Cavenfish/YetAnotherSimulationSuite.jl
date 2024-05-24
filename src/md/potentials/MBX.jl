"""
MBX potentials

This allows uing all the potentials found within the MBX package.

MBX GitHub:
https://github.com/paesanilab/MBX/tree/master

MBX Paper:
https://pubs.aip.org/aip/jcp/article/159/5/054802/2904909/MBX-A-many-body-energy-and-force-calculator-for
"""

struct _MBX_PotVars <: PotVars
  xyz
  json
  num_ats
  at_nams
  num_mon
  mon_nams
end

function MBX(bdys::Vector{Atom}, json::String)
  monomers = Dict()
  monomers["nh4+"] = ["N","H","H","H","H"]
  monomers["nh3"] = ["N","H","H","H"]
  monomers["ch4"] = ["C","H","H","H","H"]
  monomers["pf6-"] = ["P","F","F","F","F","F","F"]
  monomers["co2"] = ["C","O","O"]
  monomers["li"] = ["Li"]
  monomers["na"] = ["Na"]
  monomers["k"] = ["K"]
  monomers["rb"] = ["Rb"]
  monomers["cs"] = ["Cs"]
  monomers["f"] = ["F"]
  monomers["cl"] = ["Cl"]
  monomers["br"] = ["Br"]
  monomers["i"] = ["I"]
  monomers["ar"] = ["Ar"]
  monomers["he"] = ["He"]
  monomers["h2"] = ["H","H"]
  monomers["h2o"] = ["O","H","H"]
  monomers["n2o5"] = ["O","N","N","O","O","O","O"]

  xyz = [j for i in bdys for j in i.r]

  at_nams = [string(i.s) for i in bdys]

  mon_nams = String[]
  num_ats  = Int32[]

  tmp = values(monomers) |> (x -> length.(x)) |> unique |> sort |> reverse
  

  #TODO: add a break case to prevent infinite loop

  a = 1
  while a <= length(at_nams)
    for i in tmp
      a+i-1 > length(at_nams) && continue
      cur_ats = at_nams[a:a+i-1]
      for j in keys(monomers)
        if monomers[j] == cur_ats
          push!(mon_nams, j)
          push!(num_ats, length(monomers[j]))
          a += length(monomers[j])
          break
        end
      end
    end
  end

  num_mon = length(mon_nams)

  sym = dlsym(libmbx, :initialize_system_py_)

  _MBX_PotVars(xyz, json, num_ats, at_nams, num_mon, mon_nams)
end