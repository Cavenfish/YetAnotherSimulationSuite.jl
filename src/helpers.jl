

function CoM(bdys)
  M = sum([i.m for i in bdys])
  r = sum([i.m*i.r for i in bdys])
  return r ./ M
end


function vib_excite!(mol, eignvec, E)
  v = sqrt.((2*E) ./ mol.m) # column vector of length of atoms in mol

  mol.v .+= v .* eignvec
end
