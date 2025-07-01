function getMols(bdys::Vector{MyAtoms}, rmax; D=3)
  pts = zeros(length(bdys[1].r), length(bdys))

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

function getMols(cell::MyCell, rmax)
  n = length(cell.scaled_pos)
  p = PeriodicEuclidean(1.0)
  D = zeros(n, n)

  for i = 1:n
    for j = i+1:n
      sr = [
        p(cell.scaled_pos[i][1], cell.scaled_pos[j][1]),
        p(cell.scaled_pos[i][2], cell.scaled_pos[j][2]),
        p(cell.scaled_pos[i][3], cell.scaled_pos[j][3])
      ]
      D[i,j] = cell.lattice * sr |> norm
      D[j,i] = D[i,j]
    end
  end

  ret = dbscan(D, rmax, metric=nothing)

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