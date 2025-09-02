function scmef_getDipole(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_dipole"(pos, lat, pbc=false)
end

function scmef_getDipole(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_dipole"(pos, cell.lattice, NC=cell.NC)
end

function scmef_getInducedDipoles(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_induced_dipoles"(pos, lat, pbc=false)
end

function scmef_getInducedDipoles(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_induced_dipoles"(pos, cell.lattice, NC=cell.NC)
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

  py"scmef_get_total_electric_field"(pos, lat, pbc=false)
end

function scmef_getTotalElectricField(cell::MyCell)
  pos = getPos(cell)

  py"scmef_get_total_electric_field"(pos, cell.lattice, NC=cell.NC)
end