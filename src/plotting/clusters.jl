
function pltAlphaShape(bdys)

  pts = [i.r for i in bdys]

  A   = alphashape(pts) 

  xyz   = Tuple.(pts)
  faces = hcat(A.perimeter...) |> transpose

  # fig = Figure()
  # ax  = Axis3D(fig, xlabel=L"\AA", ylabel=L"\AA", zlable=L"\AA")

  fig = mesh(xyz, faces, color=(:purple, 0.15))


  c = [Tuple(i.r) for i in bdys if i.s=='C']
  o = [Tuple(i.r) for i in bdys if i.s=='O']
  meshscatter!(c, color=:gray, markersize=0.3)
  meshscatter!(c, color=:red,  markersize=0.3)

  fig
end 
