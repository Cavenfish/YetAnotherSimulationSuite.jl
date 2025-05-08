"""
Wrapper for phonopy python api

Notes
-------
Phonopy has strange ways of doing matrix/vector math.

To calculate the positions of atoms they do:
  [    ] |    |
   ^spos |    |
          ^lattice

Then to calculate the primitive or supercell they do:
  T * Lat

Since JMD does primitive and supercell transforms identically,
but does column vectors for spos; we cannot pass spos between
phonopy and JMD.
"""

function getNewMaskOrder(N, T)
  n::Int = det(T)
  q::Int = div(N, n)

  x = [i+(j*T[5]) for i = 1:T[1] for j = 0:T[5]-1]
  I = Int[]
  
  tmp = [i for i = 1:q:N]
  for i = 1:q
    tmp .= [j for j = i:q:N][x]
    push!(I, tmp...)
  end

  I
end


function reorderPhonopySupercell!(pos, n)
  N = length(pos)
  I = [j for i = 1:n for j = i:n:N]

  pos .= pos[I]
end

function reorderPhonopyForces!(forces, n)
  N = length(forces[1])
  I = [j for i = 1:n for j = i:n:N]
  J = [j for i = 1:n for j = i:n:N]

  for i = 1:length(I)
    J[i] = findfirst(e -> e == i, I)
  end

  for f in forces
    f .= f[J]
  end
  
end

function phonopy_getDisplacements(
  cell::MyCell, primitive, supercell; dist=0.02, symprec=1e-5)

  sFlag = sum.(supercell) |> prod

  # Imports
  phonopy      = pyimport("phonopy")
  PhonopyAtoms = phonopy.structure.atoms.PhonopyAtoms

  unitcell = PhonopyAtoms(
    symbols   = string(cell.symbols...), 
    cell      = cell.lattice, 
    positions = getPos(cell)
  )

  phonon   = phonopy.Phonopy(unitcell, 
    supercell_matrix = supercell, 
    primitive_matrix = primitive,
    factor           = phonopy.units.VaspToEv,
    symprec          = symprec
  )

  phonon.generate_displacements(distance=dist)

  supercells = phonon.supercells_with_displacements

  T        = hcat(supercell...)
  n::Int64 = det(T)
  lat      = T * cell.lattice
  mas      = repeat(cell.masses,  n)
  sym      = repeat(cell.symbols, n)
  mask     = repeat(cell.mask,    n)

  cells = MyCell[]
  for dcell in supercells
    N   = size(dcell.positions)[1]
    pos = [dcell.positions[i, :] for i = 1:N]

    if sFlag != 1
      reorderPhonopySupercell!(pos, n)
    end

    spos = getScaledPos(vcat(pos...), lat)
    vels = zero(spos)
    tmp  = Cell(lat, spos, vels, mas, sym, mask, cell.PBC, cell.NC)

    push!(cells, tmp)
  end
  
  phonon, cells
end

function phonopy_getDisplacementsDataset(
  cell::MyCell, primitive, supercell; dist=0.02, symprec=1e-5)

  # Imports
  phonopy      = pyimport("phonopy")
  PhonopyAtoms = phonopy.structure.atoms.PhonopyAtoms

  unitcell = PhonopyAtoms(
    symbols   = string(cell.symbols...), 
    cell      = cell.lattice, 
    positions = getPos(cell)
  )

  phonon   = phonopy.Phonopy(unitcell, 
    supercell_matrix = supercell, 
    primitive_matrix = primitive,
    factor           = phonopy.units.VaspToEv,
    symprec          = symprec
  )

  phonon.generate_displacements(distance=dist)

  n     = hcat(supercell...) |> det
  disps = Tuple[]
  
  for disp in phonon.displacement_dataset["first_atoms"]
    i::Int = (disp["number"] / n) + 1
    d      = disp["displacement"]
    
    push!(disps, (i,d))
  end

  phonon, disps
end

function phonopy_addForces(yamlFile, saveName, forces; symprec=1e-5)
  phonopy    = pyimport("phonopy")
  obj        = phonopy.load(yamlFile, symprec=symprec)
  n::Int64   = det(obj.supercell_matrix)

  obj.forces = forces
  obj.save(saveName)
end

function phonopy_getPhonons(phonon, path, labels, saveName; N=21)

  phonopy    = pyimport("phonopy")
  get_band_qpoints_and_path_connections = phonopy.phonon.band_structure.get_band_qpoints_and_path_connections

  phonon.produce_force_constants()
  phonon.symmetrize_force_constants()

  qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=N)

  phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
  phonon.write_yaml_band_structure(nothing, saveName)

end