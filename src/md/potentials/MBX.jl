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

function MBX(bdys)
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

  a = 1
  b = 0
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
    b += 1
    if b > length(at_nams)
      @error "Check xyz file"
    end
  end

  num_mon = length(mon_nams)

  sym = dlsym(libmbx, :initialize_system_py_)

  vars = _MBX_PotVars(xyz, MBXjson, num_ats, at_nams, num_mon, mon_nams)

  @ccall $sym(
    vars.xyz::Ptr{Cdouble},
    vars.num_ats::Ptr{Cint},
    vars.at_nams::Ptr{Ptr{Cchar}},
    vars.mon_nams::Ptr{Ptr{Cchar}},
    vars.num_mon::Ref{Cint},
    vars.json::Ptr{Cchar}
  )::Cvoid

  vars
end

function MBX(dv, v, u, p, t)

  # initialize things
  x0   = [j for i in u for j in i]
  G    = zero(x0)
  F    = Vector[]
  nats = convert(Int32, length(p.potVars.at_nams))

  E  = mbx_get_energy_grad!(G, x0, nats)

  for i = 1:3:length(G)
    push!(F, -G[i:i+2])
  end

  dv .= F ./ p.m
  if typeof(p) == NVTsimu
    p.thermostat!(p.temp,dv, v, p.m, p.thermoInps)
  end

  push!(p.energy, E)
  push!(p.forces, F)

end

function MBX(F, G, y0, p)

  nats = convert(Int32, length(p.potVars.at_nams))

  E = if G != nothing
    mbx_get_energy_grad!(G, y0, nats)
  else
    mbx_get_energy(y0, nats)
  end

  if F != nothing
    return E
  end

end

function MBX(box, scaled_pos)

  y0 = zeros(length(scaled_pos))

  for i = 1:3:length(y0)
    y0[i]   = scaled_pos[i]   * box[1]
    y0[i+1] = scaled_pos[i+1] * box[2]
    y0[i+2] = scaled_pos[i+2] * box[3]
  end

  nats = convert(Int32, length(y0)/3)
  tmp  = [box[1], 0,0,0, box[2], 0,0,0, box[3]]
  
  mbx_set_box(tmp)
  mbx_get_energy(y0, nats)
end