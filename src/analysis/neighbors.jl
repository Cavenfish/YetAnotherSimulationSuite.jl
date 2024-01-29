
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