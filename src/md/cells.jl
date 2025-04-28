#TODO:
#   -Make function to clean duplicates

struct Cell <: MyCell
  lattice::Matrix{Float64}
  scaled_pos::Vector{Vector{Float64}}
  velocity::Vector{Vector{Float64}}
  masses::Vector{Float64}
  symbols::Vector{Char}
  mask::Vector{Bool}
  PBC::Vector{Bool}
  NC::Vector{Int32}
end

function makeCell(bdys::Vector{MyAtoms}, lattice; 
  mask=repeat([false], length(bdys)), PBC=repeat([true], 3), NC=[1,1,1])

  Cell(
    lattice,
    getScaledPos(bdys, lattice),
    [i.v for i in bdys],
    [i.m for i in bdys],
    [i.s for i in bdys],
    mask, PBC, NC
  )
end

function makeBdys(cell)::Vector{MyAtoms}
  pos  = getPos(cell)
  vel  = [i for i in cell.velocity]

  [Atom(pos[i], vel[i], cell.masses[i], cell.symbols[i]) for i = 1:length(pos)]
end

function getScaledPos(x0, lattice)

  T    = inv(lattice)
  
  [T * x0[i:i+2] for i = 1:3:length(x0)-1]
end

function getPos(cell)
  [cell.lattice * i for i in cell.scaled_pos]
end

function getVolume(cell)
  a,b,c  = eachrow(cell.lattice)

  cross(b,c) |> (x -> dot(a,x))
end 

function wrap!(cell)
  r = getPos(cell)
  
  f = [floor.(transpose(cell.lattice) \ i) for i in r]

  for i = 1:length(r)
    r[i] .-= transpose(cell.lattice) * f[i]
  end

  T = inv(cell.lattice)

  for i = 1:length(r)
    cell.scaled_pos[i] .= T * r[i]
  end

end

function center!(cell)
  r  = getPos(cell)
  r0 = sum(r) ./ length(r)
  μ  = sum(0.5 * cell.lattice, dims=1) |> (x -> reshape(x, 3))
  x  = μ - r0
  
  # Translate positions
  for i in r
    i .+= x
  end
  
  # Get scaled positions
  spos = vcat(r...) |> (x -> getScaledPos(x, cell.lattice))

  # Update scaled positions
  for i in 1:length(spos)
    cell.scaled_pos[i] .= spos[i]
  end 

end

function replicate(cell, N)

  a,b,c  = eachrow(cell.lattice)
  m      = length(cell.symbols)
  n      = prod(N)
  newV   = repeat(cell.velocity, n)
  newS   = repeat(cell.symbols, n)
  newM   = repeat(cell.masses, n)
  newMa  = repeat(cell.mask, n)
  newLat = cell.lattice * Diagonal(N)
  newPos = repeat(getPos(cell), n)

  # Is this style easier to read than inline?
  f = [i*a + j*b + k*c 
        for i = 0:N[1]-1 
          for j = 0:N[2]-1 
            for k = 0:N[3]-1 
              for q = 1:m]
  
  newPos .+= f

  newScaledPos = [inv(newLat) * r for r in newPos]

  Cell(
    newLat, newScaledPos, newV, newM, newS,
    newMa, cell.PBC, cell.NC
  )
end

function replicate!(super, cell, N)

  a,b,c  = eachrow(cell.lattice)
  m      = length(cell.symbols)
  n      = prod(N)

  # Inplace update super scaled pos for non-replica atoms
  for i = 1:m
    super.scaled_pos[i] .= inv(super.lattice) * (cell.lattice * cell.scaled_pos[i])
  end

  x = m+1
  for i = 1:N[1]-1
    for j = 1:N[2]-1
      for k = 1:N[3]-1
        for q = 1:m
          # Skip if atom is masked
          mask[q] || continue

          # Get replicated pos
          r   = cell.lattice * cell.scaled_pos[q]
          r .+= (i*a) + (j*b) + (k*c)

          # inplace update super scaled pos
          super.scaled_pos[x] .= inv(super.lattice) * r
          
          # increment x
          x += 1
        end
      end
    end
  end

end

function getMIC(cell::MyCell)
  a,b,c  = eachrow(cell.lattice)
  new    = MyAtoms[]
  s      = repeat([i.s for i in bdys], 27)
  m      = repeat([i.m for i in bdys], 27)
  v      = repeat([i.v for i in bdys], 27)

  # I think it is
  f = [i*a + j*b + k*c + bdys[q].r
        for i = -1:1 
          for j = -1:1 
            for k = -1:1 
              for q = 1:length(bdys)]

  for i = 1:length(f)
    push!(new, Atom(f[i], v[i], m[i], s[i]))
  end

  makeCell(new, cell.lattice*3, mask=cell.mask, PBC=cell.PBC, NC=cell.NC)
end

function makeSuperCell(cell, T)

  lattice = T * cell.lattice

  # Clean up machine precision noise
  for i = 1:9
    if abs(lattice[i]) < 1e-8
      lattice[i] = 0.0
    end 
  end
    
  super   = replicate(cell, diag(T))
  bdys    = makeBdys(super)
  wrap!(bdys, lattice)

  makeCell(bdys, lattice, mask=cell.mask, PBC=cell.PBC, NC=cell.NC)
end

function makeSuperCell!(super, cell, T)

  super.lattice .= T * cell.lattice

  # Clean up machine precision noise
  for i = 1:9
    if abs(super.lattice[i]) < 1e-8
      super.lattice[i] = 0.0
    end 
  end
    
  replicate!(super, cell, diag(T))
  wrap!(super)

end

function getPrimitiveCell(cell, symprec)
  amu  = TOML.parsefile(joinpath(@__DIR__, "../data/Atoms.toml"))["Mass"]
  num  = TOML.parsefile(joinpath(@__DIR__, "../data/Atoms.toml"))["Number"]
  rnum = Dict(values(num) .=> keys(num))
  atms = string.(cell.symbols) |> (x -> [num[i] for i in x])
  
  L    = transpose(cell.lattice)
  tmp  = Spglib.Cell(L, cell.scaled_pos, atms)
  stan = standardize_cell(tmp, symprec)
 
  syms = [rnum[i] for i in stan.atoms]
  mas  = [ amu[i] for i in syms]
  L    = transpose(stan.lattice)

  Cell(
    L, stan.positions, cell.velocity, mas, 
    only.(syms), cell.mask, cell.PBC, cell.NC
  )
end

function getMols(cell::MyCell, rmax; D=3)
  r   = getPos(cell)

  pts = hcat(r...)

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