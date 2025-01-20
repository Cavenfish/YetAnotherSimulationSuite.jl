#TODO:
#   -Make function to center lattice
#   -Make function to get super cells
#   -Make function to get primitive cells

struct Cell
  lattice::Matrix{Float64}
  scaled_pos::Vector{Vector{Float64}}
  masses::Vector{Float64}
  symbols::Vector{Char}
end

function makeCell(bdys, lattice)

  Cell(
    lattice,
    getScaledPos(bdys, lattice),
    [i.m for i in bdys],
    [i.s for i in bdys]
  )
end

function getScaledPos(bdys, lattice)

  T    = inv(lattice)
  
  [T * i.r for i in bdys]
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

  Cell(newLat, newScaledPos, newM, newS)
end


function makeSuperCell(cell, T)

  lattice = T * cell.lattice
  n       = det(T)



end

