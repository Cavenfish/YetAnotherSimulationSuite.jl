const E_MBX2JMD    = 0.043 # kcal/mol --> eV
const MBX_GAS_JSON = joinpath(@__DIR__, "gas.json")
const MBX_PBC_JSON = joinpath(@__DIR__, "pbc.json")

# set libmbx path if MBX is properly installed
const libmbx = if haskey(ENV, "MBX_HOME")
  joinpath(ENV["MBX_HOME"], "lib/libmbx.so")
else
  @warn "MBX is not properly installed.\nMBX pot not available!"
  ""
end

function mbx_get_monomer_info(at_nams)
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

  mon_nams, num_ats
end

function mbx_get_energy(xyz, nats)

  E   = Ref{Cdouble}(0)

  ccall(
    (:get_energy_, libmbx), Cvoid,
    (Ptr{Cdouble}, Ref{Cint}, Ref{Cdouble}),
    xyz, nats, E
  )

  E[] * E_MBX2JMD
end

function mbx_get_energy_decomp(xyz, nats)
  E1b = Ref{Cdouble}(0)
  E2b = Ref{Cdouble}(0)
  E3b = Ref{Cdouble}(0)
  E4b = Ref{Cdouble}(0)
  Edp = Ref{Cdouble}(0)
  Eel = Ref{Cdouble}(0)

  ccall(
    (:get_energy_decomp_, libmbx), Cvoid,
    (
      Ptr{Cdouble}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble},
      Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}
    ),
    xyz, nats, E1b, E2b, E3b, E4b, Edp, Eel,
  )

  (E1b[], E2b[], E3b[], E4b[], Edp[], Eel[]) .* E_MBX2JMD
end

function mbx_get_energy_grad!(G, xyz, nats)

  E   = Ref{Cdouble}(0)

  ccall(
    (:get_energy_g_, libmbx), Cvoid,
    (Ptr{Cdouble}, Ref{Cint}, Ref{Cdouble}, Ptr{Cdouble}),
    xyz, nats, E, G
  )

  G .*= E_MBX2JMD

  E[] * E_MBX2JMD
end

function mbx_get_total_dipole(μ)

  ccall(
    (:get_total_dipole_, libmbx), Cvoid,
    (Ptr{Cdouble},), μ
  )

end

function mbx_get_induced_dipoles(μ)

  ccall(
    (:get_induced_dipoles_, libmbx), Cvoid,
    (Ptr{Cdouble},), μ
  )

end

function mbx_get_dipoles(μi, μp)

  ccall(
    (:get_dipoles_, libmbx), Cvoid,
    (Ptr{Cdouble}, Ptr{Cdouble},), μi, μp
  )

end

function mbx_get_charges(q)
  ccall(
    (:get_charges_, libmbx), Cvoid,
    (Ptr{Cdouble},), q
  )
end

function mbx_finalize_system()

  ccall(
    (:finalize_system_, libmbx), Cvoid,
    ()
  )

end

function mbx_set_box(box)

  l = length(box)

  ccall(
    (:set_box_, libmbx), Cvoid,
    (Ref{Cint}, Ptr{Cdouble}), l, box 
  )

end

function mbx_getDipole(cell::MyCell)
  _ = getPotEnergy(MBX(), cell)
  μ = zeros(3)

  mbx_get_total_dipole(μ)
  mbx_finalize_system()

  μ
end

function mbx_getDipole(bdys::Vector{MyAtoms})
  _ = getPotEnergy(MBX(), bdys)
  μ = zeros(3)

  mbx_get_total_dipole(μ)
  mbx_finalize_system()

  μ
end

function mbx_getDipoles(cell::MyCell)
  _  = getPotEnergy(MBX(), cell)
  N  = length(cell.masses) * 4 |> Int
  μi = zeros(N)
  μp = zeros(N)

  mbx_get_dipoles(μi, μp)
  mbx_finalize_system()

  μi = [μi[i:i+2] for i = 1:3:length(μi)]
  deleteat!(μi, 4:4:length(μi))
  μp = [μp[i:i+2] for i = 1:3:length(μp)]
  deleteat!(μp, 4:4:length(μp))

  μi, μp
end

function mbx_getDipoles(bdys::Vector{MyAtoms})
  _  = getPotEnergy(MBX(), bdys)
  N  = length(bdys) * 3 |> Int
  μi = zeros(N)
  μp = zeros(N)

  mbx_get_dipoles(μi, μp)
  mbx_finalize_system()

  μi = [μi[i:i+2] for i = 1:3:length(μi)]
  μp = [μp[i:i+2] for i = 1:3:length(μp)]

  μi, μp
end

function mbx_getConstituentEnergies(obj::Union{Vector{MyAtoms}, MyCell})
  vars = MBX(obj)
  xyz  = [j for i in vars.xyz for j in i]
  nats = convert(Int32, length(vars.at_nams))

  ret = mbx_get_energy_decomp(xyz, nats)
  mbx_finalize_system()

  ret
end

function mbx_getCharges(obj::Union{Vector{MyAtoms}, MyCell})
  vars = MBX(obj)
  q    = zeros(vars.num_mon * 4)
  
  mbx_get_charges(q)
  mbx_finalize_system()

  q
end