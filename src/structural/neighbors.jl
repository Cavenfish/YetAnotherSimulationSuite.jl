"""
    countNearestNeighbors(mol, bdys; rmax=4.5)

Count the number of neighboring molecules within a cutoff distance.

# Arguments
- `mol`: The molecule of interest.
- `bdys`: Collection of molecules to search for neighbors.
- `rmax`: (Optional) Cutoff distance for neighbor search (default 4.5).

# Returns
- Number of neighboring molecules within `rmax` of `mol`.
"""
function countNearestNeighbors(mol, bdys; rmax=4.5)
  mols = getMols(bdys)
  coms = [CoM(bdys[i]) for i in mols]
  com  = CoM(mol)
  N    = 0

  for i in coms
    norm(com - i) < rmax || continue
    N += 1
  end
  
  N
end