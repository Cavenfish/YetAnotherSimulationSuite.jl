function getMols(bdys::Vector{MyAtoms}, rmax::Float64)
  pts = zeros(length(bdys[1].r), length(bdys))

  for i = 1:length(bdys)
    pts[:, i] .= bdys[i].r
  end

  ret = dbscan(pts, rmax)

  [i.core_indices for i in ret.clusters]
end

function getPairs(bdys::Vector{MyAtoms})

  # Get mols and N
  mols = getMols(bdys, 1.5)
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

function getMols(cell::MyCell, rmax::Float64)
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
  mols = getMols(cell, 1.5)
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

function getDistanceMatrix(cell::MyCell)

  if isdiag(cell.lattice)
    D = distanceMatrixOrthorhombicCell(cell)
  else
    D = distanceMatrixAnyCell(cell)
  end

  D
end

function distanceMatrixOrthorhombicCell(cell::MyCell)
  n = length(cell.scaled_pos)
  D = zeros(n,n)
  p = PeriodicEuclidean(1.0)

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

  D
end

function distanceMatrixAnyCell(cell::MyCell)
  n     = length(cell.scaled_pos)
  D     = zeros(n,n)
  pos   = getPos(cell)
  vects = eachrow(cell.lattice)

  for i = 1:n
    for j = i+1:n
      D[i,j] = periodicDistance(pos[i], pos[j], vects)
      D[j,i] = D[i,j]
    end
  end

  D
end

function periodicDistance(x::Vector{Float64}, y::Vector{Float64}, vects)
  r = norm(x - y)

  for v in vects
    ra = norm(x - (y + v))
    rb = norm(x - (y - v))

    r  = minimum([r, ra, rb])
  end

  r
end



