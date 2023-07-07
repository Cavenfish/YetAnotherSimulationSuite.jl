export rdf, adf, density

function rdf(bdys; kwargs...)
  N = length(bdys)
  o = CoM(bdys)
  r = maximum([norm(i.r-o) for i in bdys])
  c = [CoM(bdys[i:i+1]) for i in 1:2:N]
  n = length(c)
  V = 4/3 * pi * r^3
  ρ = n / V
  d = [norm(c[i]- c[j]) for i in 1:n for j in i+1:n]
  k = kde_lscv(d; kwargs...)
  C = 4 .* pi .* ρ .* step(k.x) .* k.x .^2
  x = collect(k.x)
  y = k.density ./ C
  return x,y
end

function rdf(bdys, A; kwargs...)
  N = length(bdys)
  o = CoM(bdys)
  r = maximum([norm(i.r-o) for i in bdys])
  V = 4/3 * pi * r^3
  ρ = N / V
  a = findall(e -> e.s == A, bdys)
  n = length(a)
  d = [norm(bdys[a[i]].r - bdys[a[j]].r) for i in 1:n for j in i+1:n]
  k = kde_lscv(d; kwargs...)
  C = 4 .* pi .* ρ .* step(k.x) .* k.x .^2
  x = collect(k.x)
  y = k.density ./ C
  return x,y
end

function rdf(bdys, A, B; kwargs...)
  N = length(bdys)
  o = CoM(bdys)
  r = maximum([norm(i.r-o) for i in bdys])
  V = 4/3 * pi * r^3
  ρ = N / V
  a = findall(e -> e.s == A, bdys)
  b = findall(e -> e.s == B, bdys)
  d = [norm(bdys[i].r - bdys[j].r) for i in a for j in b]
  k = kde_lscv(d; kwargs...)
  C = 4 .* pi .* ρ .* step(k.x) .* k.x .^2
  x = collect(k.x)
  y = k.density ./ C
  return x,y
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
