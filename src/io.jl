
function readASExyz(xyz)
  sys = readlines(xyz)
  N   = length(sys) - 2
  hed = sys[2]
  sys = split.(sys[3:end], " ")
  sys = deleteat!.(sys, findall.(e -> e == "", sys))
  amu = Dict("C" => 12.011, "O" => 15.999)
  Q   = Dict("C" => -1.786, "O" => -2.332)
  set = Atom[]

  for i in range(1,N)
    props = parse.(Float64, sys[i][2:end])
    pos   = SVector{3}(props[1:3])
    s     = sys[i][1]
    q     = Q[s]

    if occursin("masses", hed)
      mas = props[4]
      vel = SVector{3}(props[5:7] ./ mas)
    else
      mas = amu[s]
      vel = SVector{3}(props[4:6] ./ mas)
    end #if-else

    particle = Atom(pos, vel, mas, q)
    push!(set, particle)
  end #for loop

  return set
end #read_ase_xyz

function write_xyz_traj(fileName::String, sr, dt)
  bodies = sr.simulation.system.bodies
  n      = length(bodies)
  i      = 0
  f      = open(fileName, "w")

  for t in sr.solution.t[1:dt:end]
    cc = get_position(sr, t)
    i += 1

    println(f, lpad(n, 9))
    println(f, "i=$i, time=$t")

    for j ∈ 1:n

      s = bodies[j].s
      x = cc[1,j]
      y = cc[2,j]
      z = cc[3,j]

      println(f, "$s   $x   $y   $z")
    end # molecule for loop

  end # time for loop
  close(f)
end #write_xyz_traj
