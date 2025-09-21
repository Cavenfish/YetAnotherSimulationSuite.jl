function scmef_getDipole(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_dipole"(pos, lat, pbc=false)
end

function scmef_getDipole(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_dipole"(pos, cell.lattice, NC=cell.NC)
end

function scmef_getTotalDipoles(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_total_dipoles"(pos, lat, pbc=false)
end

function scmef_getTotalDipoles(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_total_dipoles"(pos, cell.lattice, NC=cell.NC)
end

function scmef_getInduAndPermDipoles(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_indu_and_perm_dipoles"(pos, lat, pbc=false)
end

function scmef_getInduAndPermDipoles(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_indu_and_perm_dipoles"(pos, cell.lattice, NC=cell.NC)
end

function scmef_getConstituentEnergies(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_constituent_energies"(pos, lat, pbc=false)
end

function scmef_getConstituentEnergies(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_constituent_energies"(pos, cell.lattice, NC=cell.NC)
end

function scmef_getTotalElectricField(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_total_electric_field"(pos, lat, pbc=false) .* 27.211386
end

function scmef_getTotalElectricField(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_total_electric_field"(pos, cell.lattice, NC=cell.NC) .* 27.211386
end

function scmef_getDecomposedElectrostatics(cell::MyCell)
  factorial2(n::Int64)::Float64 = prod(n:-2:1)
  pref(rank::Int64)::Float64 = 1 / factorial2(2 * rank - 1)

  Har = 27.211386
  pos = getPos(cell)

  (
    ef, ef_d1, ef_d2, ef_d3, 
    di, di_s, qu, qu_s, oc, he
  ) = py"scmef_get_electrostatic_components"(pos, cell.lattice, NC=cell.NC)

  n = size(ef)[1]

  E_s_dip = -0.5 * pref(1) .* (di .- di_s) .* ef |> sum
  E_f_dip = 0.5 * pref(1) .* di .* ef |> sum

  E_s_quad = -0.5 * pref(2) .* (qu .- qu_s) .* ef_d1 |> sum
  E_f_quad = 0.5 * pref(2) .* qu .* ef_d1 |> sum

  E_oct = 0.5 * pref(3) .* oc .* ef_d2 |> sum
  E_hex = 0.5 * pref(4) .* he .* ef_d3 |> sum


  Dict(
    "e_self_dip" => E_s_dip * Har,
    "e_field_dip" => E_f_dip * Har,
    "e_self_quad" => E_s_quad * Har,
    "e_field_quad" => E_f_quad * Har,
    "e_field_oct" => E_oct * Har,
    "e_field_hex" => E_hex * Har
  )
end