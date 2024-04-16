

function _getDipoles4TTM!(μ, u, Q, α; μtol=1e-2, Etol=1e-6)
  μConv = false
  EConv = false
  E     = zero.(u)
  μOld  = zero.(u)
  EOld  = zero.(u)

  while !μConv || !EConv

    for i = 1:length(Q)

      #This needs to be to inplace change elements rather than
      #make pointers, or re-allocate memory
      μOld[i] .= μ[i]
      EOld[i] .= E[i]

      μ[i] .= 0.0
      E[i] .= 0.0

      for j = i+1:length(Q)

        rij = u[i] - u[j]
        r   = norm(rij)
        A   = (α[i] * α[j])^(1/6)
        c   = exp(-0.4 * (r/A)^4)
        s1  = 1 - c
        s2  = s1 - (4*0.4/3) * (r/A)^4 * c

        Eq = -s1 * rij * Q[j] / r^3
        Eμ = s2 * (3/r^5) * dot(rij, μ[j]) * rij .- s1 * μ[j] / r^3

        E[i] .+= Eq + Eμ
        μ[i] .+= α[i] * E[i]

      end      
    end

    μConv = isapprox(μOld, μ, atol=μtol) |> (x -> sum(x) == length(x))
    EConv = isapprox(EOld, E, atol=Etol) |> (x -> sum(x) == length(x))

  end

end