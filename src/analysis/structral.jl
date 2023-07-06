export rdf, adf, density

function rdf(bdys)
  
end

function adf()

end


function density(bdys, rRange)

  com = CoM(bdys)
  d   = [norm(i.r - com) for i in bdys]

  x,y = [],[]
  for r in rRange
    k = findall(e -> e <= r, d)
    M = sum([bdys[i].m for i in k]) / 6.02214e23 # convert to grams
    v = 4/3 * pi * (r * 1e-8)^3 #convert Angstrom to cm

    push!(x, r)
    push!(y, M/v)
  end
  
  return x,y
end
