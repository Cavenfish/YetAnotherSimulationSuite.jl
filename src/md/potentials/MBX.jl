"""
MBX potentials

This allows uing all the potentials found within the MBX package.

MBX GitHub:
https://github.com/paesanilab/MBX/tree/master

MBX Paper:
https://pubs.aip.org/aip/jcp/article/159/5/054802/2904909/MBX-A-many-body-energy-and-force-calculator-for
"""

MBX(;constraints=nothing) = Calculator(MBX; E=MBX, EF=MBX!, constraints=constraints)

struct _MBX_PotVars{X,J,NA,AN,NM,MN} <: PotVars
  xyz::X
  json::J
  num_ats::NA
  at_nams::AN
  num_mon::NM
  mon_nams::MN
end

function MBX(bdys::Vector{MyAtoms})

  xyz = [j for i in bdys for j in i.r]

  at_nams = [string(i.s) for i in bdys]

  mon_nams, num_ats = mbx_get_monomer_info(at_nams)

  num_mon = length(mon_nams)

  vars = if isdefined(Main, :MBX_USER_JSON)
    _MBX_PotVars(
      xyz, Main.MBX_USER_JSON, num_ats, 
      at_nams, num_mon, mon_nams
    )
  else
    _MBX_PotVars(
      xyz, MBX_GAS_JSON, num_ats, 
      at_nams, num_mon, mon_nams
    )
  end

  ccall(
    (:initialize_system_py_, libmbx), Cvoid,
    (Ptr{Cdouble}, 
     Ptr{Cint}, 
     Ptr{Ptr{Cchar}}, 
     Ptr{Ptr{Cchar}}, 
     Ref{Cint}, 
     Ptr{Cchar}
    ),
    vars.xyz, 
    vars.num_ats, 
    vars.at_nams, 
    vars.mon_nams, 
    vars.num_mon, 
    vars.json
  )

  vars
end

function MBX(cell::MyCell)

  xyz = getPos(cell)

  at_nams = [string(i) for i in cell.symbols]

  mon_nams, num_ats = mbx_get_monomer_info(at_nams)

  num_mon = length(mon_nams)

  vars = if isdefined(Main, :MBX_USER_JSON)
    _MBX_PotVars(
      xyz, Main.MBX_USER_JSON, num_ats, 
      at_nams, num_mon, mon_nams
    )
  else
    _MBX_PotVars(
      xyz, MBX_PBC_JSON, num_ats, 
      at_nams, num_mon, mon_nams
    )
  end
  
  ccall(
    (:initialize_system_py_, libmbx), Cvoid,
    (Ptr{Cdouble}, 
     Ptr{Cint}, 
     Ptr{Ptr{Cchar}}, 
     Ptr{Ptr{Cchar}}, 
     Ref{Cint}, 
     Ptr{Cchar}
    ),
    vars.xyz, 
    vars.num_ats, 
    vars.at_nams, 
    vars.mon_nams, 
    vars.num_mon, 
    vars.json
  )

  if !isdiag(cell.lattice)
    @warn "Cell lattice is not diagonal.
    MBX requires a diagonal lattice.
    Current lattice will be cast as diagonal 
      --> cell.lattice .= Diagonal(cell.lattice)"

    cell.lattice .= Diagonal(cell.lattice)
  end

  box  = reshape(cell.lattice, 9) |> (x -> convert(Vector{Float64}, x))
  mbx_set_box(box)

  vars
end

function MBX(u, p)
  x    = [j for i in u for j in i]
  nats = convert(Int32, length(p.potVars.at_nams))

  mbx_get_energy(x, nats)
end

function MBX!(F, u, p)
  x    = [j for i in u for j in i]
  G    = zero(x)
  nats = convert(Int32, length(p.potVars.at_nams))
  E    = mbx_get_energy_grad!(G, x, nats)

  for i = 1:length(F)
    j::Int = 3i - 2
    F[i]  .= -G[j:j+2]
  end
  
  E
end

function MBX(F, G, cell::MyCell, lat)
  tmp          = deepcopy(cell)
  tmp.lattice .= reshape(lat, (3,3))

  # Force all bond lengths to be 0.975 to prevent bias
  pos  = getPos(cell)
  for i = 1:3:length(pos)
    ro, rh1, rh2 = pos[[i, i+1, i+2]]
    
    rvec  = ro - rh1
    r     = norm(rvec)
    x     = 0.975 - r
    rh1 .-= x * (rvec / r)

    rvec  = ro - rh2
    r     = norm(rvec)
    x     = 0.975 - r
    rh2 .-= x * (rvec / r)

    #Maybe I can condense this somehow
  end

  # Orthogonal stress because MBX only does orthogonal
  if G != nothing
    G .= getNumericalStressOrthogonal(MBX, tmp) |> (x -> reshape(x, 9))
  end

  if F != nothing
    return getPotEnergy(MBX, tmp)
  end
end