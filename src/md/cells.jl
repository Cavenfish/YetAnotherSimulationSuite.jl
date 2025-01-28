#TODO:
#   -Make function to center lattice
#   -Make function to get primitive cells
#   -Make function to clean duplicates

struct Cell <: MyCell
  lattice::Matrix{Float64}
  scaled_pos::Vector{Vector{Float64}}
  masses::Vector{Float64}
  symbols::Vector{Char}
  PBC::Vector{Bool}
  NC::Vector{Int32}
end

function makeCell(bdys::Vector{MyAtoms}, lattice; PBC=repeat([true], 3), NC=[1,1,1])

  Cell(
    lattice,
    getScaledPos(bdys, lattice),
    [i.m for i in bdys],
    [i.s for i in bdys],
    PBC, NC
  )
end

function makeBdys(cell)
  pos  = getPos(cell)
  vel  = [zeros(3) for i in pos]

  [Atom(pos[i], vel[i], cell.masses[i], cell.symbols[i]) for i = 1:length(pos)]
end

function getScaledPos(bdys::Vector{MyAtoms}, lattice)

  T    = inv(lattice)
  
  [T * i.r for i in bdys]
end

function getScaledPos(x0, lattice)

  T    = inv(lattice)
  
  [T * x0[i:i+2] for i = 1:3:length(x0)-1]
end

function getPos(cell)
  [cell.lattice * i for i in cell.scaled_pos]
end

function wrap!(bdys, lattice)

  f = [floor.(transpose(lattice) \ i.r) for i in bdys]
  
  for i = 1:length(bdys)
    bdys[i].r .-= transpose(lattice) * f[i]
  end

end

# function center!()
# end

function replicate(cell, N)

  a,b,c  = eachrow(cell.lattice)
  m      = length(cell.symbols)
  n      = prod(N)
  newS   = repeat(cell.symbols, n)
  newM   = repeat(cell.masses, n)
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

  Cell(newLat, newScaledPos, newM, newS, cell.PBC)
end


function makeSuperCell(cell, T)

  lattice = T * cell.lattice
  N       = [T[1,1], T[2,2], T[3,3]]
  
  super   = replicate(cell, N)
  bdys    = makeBdys(super)
  wrap!(bdys, lattice)

  makeCell(bdys, lattice, PBC=cell.PBC)
end

