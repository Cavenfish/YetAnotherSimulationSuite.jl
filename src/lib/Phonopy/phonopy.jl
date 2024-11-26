
struct PhonopyCellInfo
  cell::Matrix{Float64}
  scaled_positions::Vector{Vector{Float64}}
  supercell_matrix::Vector{Vector{Float64}}
  primitive_matrix::Vector{Vector{Float64}}
end

function phonopy_generate_displacements(syms, cellInfo, dist)
  phonopy  = pyimport("phonopy")
  Atoms    = pyimport("phonopy.structure.atoms")

  unitcell = Atoms.PhonopyAtoms(symbols=syms, cell=cellInfo.cell, scaled_positions=cellInfo.scaled_positions)
  phonon   = phonopy.Phonopy(unitcell, supercell_matrix=cellInfo.supercell_matrix, primitive_matrix=cellInfo.primitive_matrix)

  phonon.generate_displacements(distance=dist)                          

  phonon.supercells_with_displacements
end

function phonopy_save(obj, forces)
  obj.forces = forces
  obj.produce_force_constants()
  obj.save()
end

function phonopy_load_and_get_DOS(yamlFile, path, labels, npoints)
  phonopy  = pyimport("phonopy")

  phonon = phonopy.load(yamlFile)

  qpoints, connections = phonopy.phonon.band_structure.get_band_qpoints_and_path_connections(path, npoints=npoints)


  phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)

  # To plot DOS next to band structure
  phonon.run_mesh([20, 20, 20])
  phonon.run_total_dos()
  phonon.write_total_dos()

  # To plot PDOS next to band structure
  phonon.run_mesh([20, 20, 20], with_eigenvectors=True, is_mesh_symmetry=False)
  phonon.run_projected_dos()
  phonon.write_projected_dos()

end