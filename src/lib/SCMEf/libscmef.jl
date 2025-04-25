function scmef_get_dipole(bdys::Vector{MyAtoms})
  pos = [i.r for i in bdys]
  lat = [100 0 0;0 100 0; 0 0 100]

  py"scmef_get_dipole"(pos, lat, pbc=false)
end

function scmef_get_dipole(cell::MyCell)

  pos = getPos(cell)

  py"scmef_get_dipole"(pos, cell.lattice, NC=cell.NC)
end