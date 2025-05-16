"""
JMD file IO is done using Chemfiles

Source repo:
https://github.com/chemfiles/chemfiles/tree/master

Julia bindings repo:
https://github.com/chemfiles/Chemfiles.jl/tree/master
"""

function readSystem(file::String)
  buf    = Trajectory(file)
  N::Int = length(buf)

  if N > 1
    return readTraj(buf, N)
  else
    return read(buf) |> readFrame
  end

end

function readFrame(frame::Frame)
  N::Int = length(frame)
  bdys   = MyAtoms[]
  pos    = positions(frame)
  n      = size(pos)[1]
  ucell  = UnitCell(frame)
  lat    = matrix(ucell)
  
  has_velocities(frame) ? vel = velocities(frame) : vel = zero(pos)

  atm   = Atom(frame, 0)
  props = list_properties(atm)

  for i = 1:N
    # C++ Indexing (ie. starts at 0)
    atm = Atom(frame, i-1)

    if "velo" in props
      vel[:, i] .= property(atm, "velo")
    end
    
    r = MVector{n}(pos[:, i])
    v = MVector{n}(vel[:, i])
    m = mass(atm)
    s = name(atm)

    push!(bdys, Particle(r, v, m, s))
  end

  sum(lat) == 0.0 && return bdys

  makeCell(bdys, lat)
end

function readTraj(buf::Trajectory, N::Int)
  zeros(N)
end

function Base.write(file::String, bdys::Vector{MyAtoms})
  buf   = Trajectory(file, 'w')
  frame = Frame()
  r     = zeros(3)
  v     = zeros(3)

  add_velocities!(frame)

  for i in bdys
    r  .= i.r
    v  .= i.v
    atm = Atom(i.s)
    
    if occursin(".xyz", file)
      set_property!(atm, "velo", v)
      add_atom!(frame, atm, r)
    else
      add_atom!(frame, atm, r, v)
    end
  end

  write(buf, frame)

  close(buf)
end

function Base.write(file::String, cell::MyCell)
  buf   = Trajectory(file, 'w')
  frame = Frame()
  ucell = UnitCell(cell.lattice)
  pos   = getPos(cell)

  set_cell!(buf, ucell)
  add_velocities!(frame)
  
  for i = 1:length(cell.masses)
    r = pos[i]
    v = cell.velocity[i]
    s = cell.symbols[i]

    add_atom!(frame, Atom(s), r, v)
  end

  write(buf, frame)

  close(buf)
end