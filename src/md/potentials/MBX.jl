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

function MBX(bdys::Vector{MyAtoms})

  xyz = [j for i in bdys for j in i.r]

  at_nams = [string(i.s) for i in bdys]

  mon_nams, num_ats = mbx_get_monomer_info(at_nams)

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

function MBX(cell::MyCell)

  xyz = getPos(cell)

  at_nams = [string(i) for i in cell.symbols]

  mon_nams, num_ats = mbx_get_monomer_info(at_nams)

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

  box  = reshape(cell.lattice, 9) |> (x -> convert(Vector{Float64}, x))
  mbx_set_box(box)

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

# MBX can only handle Diagonal lattices
function MBX(L, scaled_pos, p)

  println(L)
  # L    = reshape(box, (3,3)) |> Diagonal
  r    = [L * i for i in scaled_pos]
  y0   = [j for i in r for j in i]
  nats = length(r)
  box  = reshape(L, 9) |> (x -> convert(Vector{Float64}, x))
  
  mbx_set_box(box)
  mbx_get_energy(y0, nats)
end