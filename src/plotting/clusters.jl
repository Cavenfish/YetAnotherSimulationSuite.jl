
function pltAlphaShape(bdys)

  pts = [i.r for i in bdys]

  A   = alphashape(pts) 

  xyz   = Tuple.(pts)
  faces = hcat(A.perimeter...) |> transpose

  fig = Figure()

  mesh(xyz, faces, color=(:purple, 0.15))

  fig
end 
