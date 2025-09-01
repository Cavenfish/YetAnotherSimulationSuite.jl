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

"""
    ifProperty(frame, props, p)

Return the value of property `p` from `frame` if it exists in `props`, else return 0.0.

# Arguments
- `frame`: Chemfiles Frame object.
- `props`: List of property names.
- `p`: Property name to look for.

# Returns
- Value of the property as Float64, or 0.0 if not found.
"""
function ifProperty(frame, props, p)
  if p in props
    x = property(frame, p)
    return parse(Float64, x)
  else 
    return 0.0
  end
end

"""
    atomForces(frame, i)

Get the force vector for atom `i` in the given `frame`.

# Arguments
- `frame`: Chemfiles Frame object.
- `i`: Atom index.

# Returns
- Static vector of forces.
"""
function atomForces(frame, i)
  x = property(Atom(frame, i), "forces")
  n = length(x)

  SVector{n}(x)
end

"""
    readSystem(file::String)

Read a system from a file and return either a trajectory or a single frame.

# Arguments
- `file`: Path to the file.

# Returns
- Traj object if multiple frames, otherwise a single frame as MyCell or MyAtoms.
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

"""
    readFrame(frame::Frame)

Convert a Chemfiles Frame to a vector of MyAtoms or a MyCell if lattice is present.

# Arguments
- `frame`: Chemfiles Frame object.

# Returns
- Vector of MyAtoms or MyCell object.
"""
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

"""
    readImage(frame::Frame)

Read an image (snapshot) from a Chemfiles Frame, extracting positions, velocities, forces, and properties.

# Arguments
- `frame`: Chemfiles Frame object.

# Returns
- Image object containing positions, velocities, time, temperature, energy, and forces.
"""
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

"""
    readTraj(buf::Trajectory, N::Int)

Read a trajectory from a Chemfiles Trajectory buffer.

# Arguments
- `buf`: Chemfiles Trajectory object.
- `N`: Number of frames to read.

# Returns
- Traj object containing all images, masses, symbols, and lattice.
"""
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

"""
    Base.write(file::String, bdys::Vector{MyAtoms})

Write a vector of MyAtoms to a file using Chemfiles.

# Arguments
- `file`: Output file path.
- `bdys`: Vector of MyAtoms to write.

# Side Effects
- Writes atomic coordinates and velocities to file.
"""
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

"""
    Base.write(file::String, cell::MyCell)

Write a MyCell object to a file using Chemfiles.

# Arguments
- `file`: Output file path.
- `cell`: MyCell object to write.

# Side Effects
- Writes atomic coordinates, velocities, and cell information to file.
"""
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

"""
    Base.write(file::String, traj::MyTraj; step=1)

Write a trajectory to a file using Chemfiles.

# Arguments
- `file`: Output file path.
- `traj`: MyTraj object to write.
- `step`: Step interval for writing frames (default: 1).

# Side Effects
- Writes trajectory frames, including positions, velocities, energies, and forces, to file.
"""
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

    # Kinectic Energy
    K = 0.5 * sum( traj.masses .* dot.(img.vel, img.vel))

    # Add frame properties
    set_property!(frame, "time", img.time)
    set_property!(frame, "temp", img.temp)
    set_property!(frame, "energy", img.energy + K)

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