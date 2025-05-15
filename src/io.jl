
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

  for i = 1:N
    # C++ Indexing (ie. starts at 0)
    atm = Atom(frame, i-1)
    
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
    r .= i.r
    v .= i.v

    add_atom!(frame, Atom(i.s), r, v)
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

function readCell(fileName::String)
  stream = readlines(fileName)
  amu    = TOML.parsefile(joinpath(@__DIR__, "data/Atoms.toml"))["Mass"]
  bdys   = MyAtoms[]

  tmp  = split(stream[2], "Lattice=")[2] |> (
    x -> split(x, "\"")[2]) |> (x -> split(x, " "))
  lat  = parse.(Float64, tmp) |> (x -> reshape(x, (3,3))) |> transpose
  
  for line in stream[3:end]
    s = split(line, " ")
    s = deleteat!(s, findall(e -> e == "", s))

    pos  = parse.(Float64, s[2:4])
    pos  = Vector(pos)
    
    vel  = if length(s) >= 7
      tmp = parse.(Float64, s[5:end])
      Vector(tmp)
    else
      zeros(3)
    end
    
    mas  = amu[s[1]]
    sym  = s[1][1]

    atom = Atom(pos, vel, mas, sym)
    push!(bdys, atom)
  end

  makeCell(bdys, lat)
end

function writeCell(fileName::String, cell)
  f = open(fileName, "w")
  N = length(cell.symbols)
  
  l = replace("$(cell.lattice)", [';', '[', ']'] => "")
  r = getPos(cell)

  println(f, N)
  println(f, "Lattice=\"$l\"")

  for i in 1:N

    s          = cell.symbols[i]
    x,y,z      = r[i]

    println(f, "$s   $x   $y   $z")
  end 

  close(f)
end

function writeXyzTraj(fileName::String, solu; dt=1, lat=nothing)
  f    = open(fileName, "w")
  bdys = solu.prob.p.bdys
  N    = length(bdys)
  T    = length(solu.t)

  for i in 1:dt:T
    t = solu.t[i]
    u = solu.u[i].x[2] # x[1] -> vel || x[2] -> pos

    println(f, N)

    if lat != nothing
      l = replace("$(lat)", [';', '[', ']'] => "")
      println(f, "i=$i, time=$t, Lattice=\"$l\"")
    else
      println(f, "i=$i, time=$t")
    end

    for j in 1:N

      s          = bdys[j].s
      x,y,z      = u[j]
      vx, vy, vz = bdys[j].v 

      println(f, "$s   $x   $y   $z   $vx   $vy   $vz")
    end 

  end
  close(f)
end 

function writeXyzTraj(fileName::String, tj::MyTraj; dt=1, lat=nothing)
  f = open(fileName, "w")
  T = length(tj.t)
  N = length(tj.m)

  for i in 1:dt:T
    t = tj.t[i]
    u = tj.r[i]

    println(f, N)

    if lat != nothing
      l = replace("$(lat)", [';', '[', ']'] => "")
      println(f, "i=$i, time=$t, Lattice=\"$l\"")
    else
      println(f, "i=$i, time=$t")
    end

    for j in 1:N

      s     = tj.s[j]
      x,y,z = u[j]

      println(f, "$s   $x   $y   $z")
    end 

  end
  close(f)
end  

