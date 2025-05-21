function getMols(bdys::Vector{MyAtoms}, rmax; D=3)
  pts = zeros(length(r[1]), length(r))

  for i = 1:length(bdys)
    pts[:, i] .= bdys[i].r
  end

  ret = dbscan(pts[1:D, :], rmax)

  [i.core_indices for i in ret.clusters]
end

function getPairs(bdys::Vector{MyAtoms})

  # Get mols and N
  mols = if length(bdys) <= 3
    getMols(bdys, 1.5, D=length(bdys)-1) 
  else
    getMols(bdys, 1.5)
  end
  N    = size(mols)[1]

  # Make all pairs
  pars = Pair[]
  for i in 1:N
    for j in i+1:N
      push!(pars, Pair(mols[i],mols[j]))
    end
  end

  pars, mols
end

function getMols(cell::MyCell, rmax; D=3)
  r   = getPos(cell)
  pts = zeros(length(r[1]), length(r))

  for i = 1:length(r)
    pts[:, i] .= r[i]
  end

  ret = dbscan(pts[1:D, :], rmax)

  [i.core_indices for i in ret.clusters]
end

function getPairs(cell::MyCell)

  n = length(cell.masses)

  # Get mols and N
  mols = if n <= 3
    getMols(cell, 1.5, D=n-1) 
  else
    getMols(cell, 1.5)
  end
  N    = size(mols)[1]

  # Make all pairs
  pars = Pair[]
  for i in 1:N
    for j in i+1:N
      push!(pars, Pair(mols[i],mols[j]))
    end
  end

  pars, mols
end