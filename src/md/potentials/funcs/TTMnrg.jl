

function _interTTM!(F, u, μ, i, j, Qi, Qj, Aij, Bij, Cij, A6ij; kwargs...)

  E  = 0.0
  E += _shortDisp!(F, u, i, j, Aij, Bij)
  E += _longDisp!( F, u, i, j, Cij; kwargs...)
  E += _Vpol4Fcc!( F, u, i, j, Qi, Qj, A6ij)
  _Vpol4Fcd!( F, u, i, j, Qi, Qj, μ[i], μ[j], A6ij)
  _Vpol4Fdd!( F, u, i, j, Qi, Qj, μ[i], μ[j], A6ij)
  _Vpol4Fcd!( F, u, j, i, Qj, Qi, μ[j], μ[i], A6ij)
  _Vpol4Fdd!( F, u, j, i, Qj, Qi, μ[j], μ[i], A6ij)

  E
end

function _getDipolePolarizationEnergy(μ, α)
  E = 0.0
  for i = 1:length(α)
    E += dot(μ[i], μ[i]) / (2 * α[i])
  end
  E
end

function _getPermanentEfield!(Eq, u, Q, α)

  for i = 1:length(u)
    Eq[i] .*= 0.0
    for j = 1:length(u)

      i == j && continue

      rij = u[i] - u[j]
      r   = norm(rij)
      A   = (α[i] * α[j])^(1/6)
      c   = exp(-0.4 * (r/A)^4)
      s1  = 1 - c

      Eq[i] .+= -s1 * rij * Q[j] / r^3
    end
  end

end

function _getDipoles4TTM_MatrixInversion!(μ, u, Q, α, Eq, M)

  Tij(r, A, a) = (1 - exp(-a*(r/A)^4) + (a^(1/4) * r / A) * gamma(3/4, a * (r/A)^4)) / r

  n  = length(α)

  for i = 1:n
    for j = 1:n

      if i == j
        M[i,j] = α[i]^-1
      else
        r      = u[i] - u[j] |> norm
        A      = (α[i]*α[j])^(1/6)
        M[i,j] = - Tij(r, A, 0.055)
      end

    end
  end

  tmp = inv(M) * Eq

  #Inplace swap mu values
  for i = 1:length(tmp)
    μ[i] .= tmp[i]
  end

end

function _Vpol4Fcc(ri, rj, Qi, Qj, A)
  rvec = rj - ri
  r    = norm(rvec)
  a    = 0.4

  c    = a * (r/A)^4
  s0   = 1 - exp(-c) + (a^0.25 * r / A) * gamma(0.75, c)
  s1   = 1 - exp(-c)

  Eqq  = s0 * Qi * Qj / r
  Fqq  = s1 * Qi * Qj / r^3 * rvec

  Eqq, Fqq
end

function _Vpol4Fcc!(F, u, i, j, Qi, Qj, A)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  a     = 0.4

  c     = a * (r/A)^4
  s0    = 1 - exp(-c) + (a^0.25 * r / A) * gamma(0.75, c)
  s1    = 1 - exp(-c)

  Eqq   = s0 * Qi * Qj / r
  Fqq   = s1 * Qi * Qj / r^3 * rvec

  F[i] -= Fqq
  F[j] += Fqq

  Eqq
end

function _Vpol4Fcd(ri, rj, Qi, Qj, μi, μj, A)
  rvec = rj - ri
  r    = norm(rvec)
  a    = 0.4

  c     = exp(-a*(r/A)^4)
  s1    =  1 - c
  s2    = s1 - (4a/3) * (r/A)^4 * c
  s3    = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi   = dot(rvec, μi)
  rμj   = dot(rvec, μj)
  Equ   = (Qi * rμj - Qj * rμi) * (s1 / r^3)
  Fqu   = (Qi * rμj - Qj * rμi) * (s2 * 3 / r^5 * rvec) .+ (Qj * μi .- Qi * μj) * (s1 / r^3)

  Equ, Fqu
end

function _Vpol4Fcd!(F, u, i, j, Qi, Qj, μi, μj, A)
  rvec  = u[j] - u[i]
  r     = norm(rvec)
  a     = 0.4

  c     = exp(-a*(r/A)^4)
  s1    =  1 - c
  s2    = s1 - (4a/3) * (r/A)^4 * c
  s3    = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi   = dot(rvec, μi)
  rμj   = dot(rvec, μj)
  Fqu   = (Qi * rμj - Qj * rμi) * (s2 * 3 / r^5 * rvec) .+ (Qj * μi .- Qi * μj) * (s1 / r^3)

  F[i] -= Fqu

end

function _Vpol4Fdd(ri, rj, Qi, Qj, μi, μj, A)
  rvec = rj - ri
  r    = norm(rvec)
  a    = 0.055

  c    = exp(-a*(r/A)^4)
  s1   =  1 - c
  s2   = s1 - (4a/3) * (r/A)^4 * c
  s3   = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi  = dot(rvec, μi)
  rμj  = dot(rvec, μj)
  μij  = dot(μi, μj)

  Euu  = (s1 / r^3 * μij) - (s2 * 3 / r^5 * rμi * rμj)
  Fuu  = (-s3 * 15 / r^7 * rμi * rμj * rvec) .+ (s2 * 3 / r^5) * (μij * rvec .+ rμi * μj .+ rμj * μi)

  Euu, Fuu
end

function _Vpol4Fdd!(F, u, i, j, Qi, Qj, μi, μj, A)
  rvec = u[j] - u[i]
  r    = norm(rvec)
  a    = 0.055

  c    = exp(-a*(r/A)^4)
  s1   =  1 - c
  s2   = s1 - (4a/3) * (r/A)^4 * c
  s3   = s2 - (4a/15) * (r/A)^4 * (4a * (r/4)^4 - 1) * c

  rμi  = dot(rvec, μi)
  rμj  = dot(rvec, μj)
  μij  = dot(μi, μj)

  Fuu  = (-s3 * 15 / r^7 * rμi * rμj * rvec) .+ (s2 * 3 / r^5) * (μij * rvec .+ rμi * μj .+ rμj * μi)

  F[i] -= Fuu
  
end


#TODO:
  # Make iterative dipole function have a cutoff

  function _getDipoles4TTM_Iterative!(μ, u, Q, α, mols, Eq; μtol=1e-2, Etol=1e-6)
    μConv = false
    EConv = false
    E     = zero.(μ)
    μOld  = zero.(μ)
    EOld  = zero.(μ)
  
    iter = 1
    while !μConv || !EConv
  
      for i = 1:length(μ)
        #This needs to be here to have inplace change elements
        #rather than make pointers, or re-allocate memory
        μOld[i] .= μ[i]
        EOld[i] .= E[i]
  
        μ[i] .= 0.0
        E[i] .= 0.0
      end
  
      # iter > 5 && break
  
      for n = 1:length(mols)
  
        for m = 1:length(mols)
  
          n == m && continue
  
          for i in mols[n]
            for j in mols[m]
  
              rij = u[i] - u[j]
              r   = norm(rij)
              A   = (α[i] * α[j])^(1/6)
              c   = exp(-0.055 * (r/A)^4)
              s1  = 1 - c
              s2  = s1 - (4*0.055/3) * (r/A)^4 * c
      
              Eμ = s2 * (3/r^5) * dot(rij, μOld[j]) * rij .- s1 * μOld[j] / r^3
      
              E[i] .+= Eq[i] + Eμ
              μ[i] .+= α[i] * E[i]
  
            end
          end
        end
      end
  
      # println(μOld[1])
      # println(μ[1])
      # (μOld .- μ) |> maximum |> println
      # (μOld .- μ) |> (x -> findfirst(e -> e == maximum(x), x)) |> println
  
      μConv = isapprox(μOld, μ, atol=μtol) |> (x -> sum(x) == length(x))
      EConv = isapprox(EOld, E, atol=Etol) |> (x -> sum(x) == length(x))
      iter +=1
    end
  
    println("Total Iterations $iter")
  
  end

  