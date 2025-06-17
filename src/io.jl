"""
JMD file IO is done using Chemfiles

Source repo:
https://github.com/chemfiles/chemfiles/tree/master

Julia bindings repo:
https://github.com/chemfiles/Chemfiles.jl/tree/master
"""

atomMass(frame::Frame, i::Int) = Atom(frame, i) |> mass
atomName(frame::Frame, i::Int) = Atom(frame, i) |> name

getMasses(frame::Frame) = [atomMass(frame, i) for i = 0:length(frame)-1]
getNames( frame::Frame) = [atomName(frame, i) for i = 0:length(frame)-1]

function ifProperty(frame, props, p)
  if p in props
    x = property(frame, p)
    return parse(Float64, x)
  else 
    return 0.0
  end
end

function atomForces(frame, i)
  x = property(Atom(frame, i), "forces")
  n = length(x)

  SVector{n}(x)
end

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

function readImage(frame::Frame)
  pos    = positions(frame)
  props  = list_properties(frame)
  n      = size(pos)[1]

  forces = if "forces" in list_properties(Atom(frame, 1))
    [atomForces(frame, i) for i = 0:length(frame)-1]
  else
    tmp = zero(pos)
    [SVector{n}(tmp[:, i]) for i = 1:length(frame)]
  end

  # Props here must be able to be Float
  time   = ifProperty(frame, props, "time")
  temp   = ifProperty(frame, props, "temp")
  energy = ifProperty(frame, props, "energy")

  # Check for velocities
  has_velocities(frame) ? vel = velocities(frame) : vel = zero(pos)

  # Reshape pos & vel
  pos = [SVector{n}(pos[:, i]) for i = 1:length(frame)]
  vel = [SVector{n}(vel[:, i]) for i = 1:length(frame)]

  Image(pos, vel, time, temp, energy, forces)
end

function readTraj(buf::Trajectory, N::Int)
  frame   = read(buf)
  lattice = UnitCell(frame) |> matrix
  symbols = getNames(frame)
  masses  = getMasses(frame)
  images  = [readImage(frame)]

  for i = 2:N
    frame = read(buf)
    push!(images, readImage(frame))
  end

  Traj(images, masses, symbols, lattice)
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
  ucell = Matrix(cell.lattice) |> UnitCell
  pos   = getPos(cell)
  r     = zeros(3)
  v     = zeros(3)

  set_cell!(buf, ucell)
  add_velocities!(frame)
  
  for i = 1:length(cell.masses)
    r  .= pos[i]
    v  .= cell.velocity[i]
    s   = cell.symbols[i]
    atm = Atom(s)

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

function Base.write(file::String, traj::MyTraj; step=1)
  buf   = Trajectory(file, 'w')
  ucell = Matrix(traj.lattice) |> UnitCell
  r     = zeros(3)
  v     = zeros(3)
  f     = zeros(3)

  set_cell!(buf, ucell)

  for img in traj.images[1:step:end]
    frame = Frame()
    add_velocities!(frame)

    # Add frame properties
    set_property!(frame, "time", img.time)
    set_property!(frame, "temp", img.temp)
    set_property!(frame, "energy", img.energy)

    for i = 1:length(img.pos)
      r  .= img.pos[i]
      v  .= img.vel[i]
      f  .= img.forces[i]
      atm = Atom(traj.symbols[i])

      set_property!(atm, "forces", f)
      
      if occursin(".xyz", file)
        set_property!(atm, "velo", v)
        add_atom!(frame, atm, r)
      else
        add_atom!(frame, atm, r, v)
      end
    end

    write(buf, frame)
  end

  close(buf)
end