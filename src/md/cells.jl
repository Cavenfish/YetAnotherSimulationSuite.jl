#TODO:
#   -Make function to clean duplicates

struct Cell{D, B, I<:Int, F<:AbstractFloat, S<:AbstractString} <: MyCell
  lattice::MMatrix{D,D,F}
  scaled_pos::Vector{MVector{D,F}}
  velocity::Vector{MVector{D,F}}
  masses::Vector{F}
  symbols::Vector{S}
  mask::Vector{B}
  PBC::Vector{B}
  NC::Vector{I}
end

function Cell(
  lat::Matrix, spos::AbstractArray, vel::AbstractArray, 
  mas::Vector{F}, sym::Vector{S}, mask::Vector{Bool},
  PBC::Vector{Bool}, NC::Vector{Int}
) where {F<:AbstractFloat, S<:AbstractString}

  n = length(spos[1])

  Cell(
    MMatrix{n,n}(lat),
    MVector{n}.(spos),
    MVector{n}.(vel),
    mas, sym, mask, PBC, NC
  )
end

function makeCell(bdys::Vector{MyAtoms}, lattice::AbstractMatrix; 
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

function trim!(cell::MyCell, iter)
  deleteat!(cell.scaled_pos, iter)
  deleteat!(cell.velocity, iter)
  deleteat!(cell.masses, iter)
  deleteat!(cell.symbols, iter)
  deleteat!(cell.mask, iter)
end

function makeBdys(cell::MyCell)::Vector{MyAtoms}
  D    = length(cell.velocity[1])
  pos  = MVector{D}.(getPos(cell))
  vel  = [MVector{D}(i) for i in cell.velocity]

  [Particle(pos[i], vel[i], cell.masses[i], cell.symbols[i]) for i = 1:length(pos)]
end

function getScaledPos(x0, lattice)
  T = inv(lattice)
  
  [T * x0[i:i+2] for i = 1:3:length(x0)-1]
end

function getScaledPos!(cell, x0)
  T = inv(cell.lattice)

  for i = 1:3:length(x0)-1
    j::Int = (i+2)/3

    cell.scaled_pos[j] .= T * x0[i:i+2]
  end

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

function Base.repeat(cell::MyCell, count::Integer)
  Cell(
    deepcopy(cell.lattice),
    myRepeat(cell.scaled_pos, count, cell.mask),
    myRepeat(cell.velocity, count, cell.mask),
    myRepeat(cell.masses, count, cell.mask),
    myRepeat(cell.symbols, count, cell.mask),
    myRepeat(cell.mask, count, cell.mask),
    deepcopy(cell.PBC), 
    deepcopy(cell.NC)
  )
end

function replicate!(super, cell, N)

  a,b,c  = eachrow(cell.lattice)
  m      = length(cell.symbols)
  M      = length(super.scaled_pos)
  n      = prod(N) * m
  l      = n - (sum(cell.mask) * (prod(N) - 1))

  # End if super is too small 
  # TODO: make this increase super size
  if M < l 
    @error "super is too small"
    return
  end

  # Trim super (if needed)
  if M > l
    trim!(super, l+1:M)
  end

  x = 1
  for i = 0:N[1]-1
    for j = 0:N[2]-1
      for k = 0:N[3]-1
        for q = 1:m
          # Skip replication if atom is masked
          cell.mask[q] && x > m && continue

          # Get replicated pos
          r   = cell.lattice * cell.scaled_pos[q]
          r .+= (i*a) + (j*b) + (k*c)

          # inplace update super
          super.scaled_pos[x] .= inv(super.lattice) * r
          super.velocity[x]   .= cell.velocity[q]
          super.masses[x]      = cell.masses[q]
          super.symbols[x]     = cell.symbols[q]
          super.mask[x]        = cell.mask[q]
          
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
    push!(new, Particle(f[i], v[i], m[i], s[i]))
  end

  makeCell(new, cell.lattice*3, mask=cell.mask, PBC=cell.PBC, NC=cell.NC)
end

function makeSuperCell(cell, T)

  N       = diag(T)
  lattice = T * cell.lattice

  # Clean up machine precision noise
  for i = 1:9
    if abs(lattice[i]) < 1e-8
      lattice[i] = 0.0
    end 
  end

  # allocate super cell
  super = repeat(cell, prod(N))

  # inplace update lattice
  super.lattice .= lattice

  # inplace update fields
  replicate!(super, cell, N)

  # inplace wrap atoms outside PBC
  wrap!(super)

  super
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