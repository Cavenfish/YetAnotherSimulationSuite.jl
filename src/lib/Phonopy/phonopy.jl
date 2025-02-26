"""
Wrapper for phonopy python api
"""

# For now supercell is fixed, can become variable if the 
# reordering problem is solved
function phonopy_getDisplacements(cell::MyCell, saveName, primitive; 
                                  dist=0.02, symprec=1e-5)
  # fixed supercell matrix
  supercell = [[1, 0, 0], [0, 1, 0,], [0, 0, 1]]

  # Imports
  phonopy      = pyimport("phonopy")
  PhonopyAtoms = phonopy.structure.atoms.PhonopyAtoms


  unitcell = PhonopyAtoms(
    symbols          = string(cell.symbols...), 
    cell             = cell.lattice, 
    scaled_positions = cell.scaled_pos
  )
  phonon   = phonopy.Phonopy(unitcell, 
    supercell_matrix = supercell, 
    primitive_matrix = primitive,
    factor           = phonopy.units.VaspToEv,
    symprec          = symprec
  )

  phonon.generate_displacements(distance=dist)

  phonon.save(saveName)

  supercells = phonon.supercells_with_displacements

  cells = MyCell[]
  for dcell in supercells
    tmp = deepcopy(cell)

    for i = 1:length(tmp.scaled_pos)
      tmp.scaled_pos[i] .= dcell.scaled_positions[i, :]
    end

    push!(cells, tmp)
  end
  
  phonon, cells
end

function phonopy_addForces(yamlFile, saveName, forces; symprec=1e-5)
  phonopy    = pyimport("phonopy")
  obj        = phonopy.load(yamlFile, symprec=symprec)
  obj.forces = forces
  obj.save(saveName)
end

function phonopy_getPhonons(phonon, path, labels; N=21)

  phonopy    = pyimport("phonopy")
  get_band_qpoints_and_path_connections = phonopy.phonon.band_structure.get_band_qpoints_and_path_connections

  phonon.produce_force_constants()
  phonon.symmetrize_force_constants()
  phonon.save("params.yaml")

  qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=N)

  phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
  phonon.write_yaml_band_structure(None, "band.yaml")

end