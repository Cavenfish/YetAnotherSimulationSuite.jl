
#TODO:
  # Make iterative dipole function have a cutoff


function _interTTM!(F, u, μ, i, j, Qi, Qj, Aij, Bij, Cij, A6ij; kwargs...)

  E  = 0.0
  E += _Coulomb!(  F, u, i, j, Qi, Qj)
  E += _shortDisp!(F, u, i, j, Aij, Bij)
  E += _longDisp!( F, u, i, j, Cij; kwargs...)
  E += _Vpol4Fcc!( F, u, i, j, Qi, Qj, A6ij)
  E += _Vpol4Fcd!( F, u, i, j, Qi, Qj, μ[i], μ[j], A6ij)
  E += _Vpol4Fdd!( F, u, i, j, Qi, Qj, μ[i], μ[j], A6ij)

  E += _Vpol4Fcc!( F, u, j, i, Qj, Qi, A6ij)
  E += _Vpol4Fcd!( F, u, j, i, Qj, Qi, μ[j], μ[i], A6ij)
  E += _Vpol4Fdd!( F, u, j, i, Qj, Qi, μ[j], μ[i], A6ij)

  E
end

function _getDipoles4TTM!(μ, u, Q, α, mols; μtol=1e-6, Etol=1e-6)
  μConv = false
  EConv = false
  E     = zero.(μ)
  μOld  = zero.(μ)
  EOld  = zero.(μ)

  while !μConv || !EConv

    for n = 1:length(mols)

      for m = n+1:length(mols)

        for i in mols[n]
          #This needs to be to inplace change elements rather than
          #make pointers, or re-allocate memory
          μOld[i] .= μ[i]
          EOld[i] .= E[i]
    
          μ[i] .= 0.0
          E[i] .= 0.0
          for j in mols[m]

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
      end
    end

    μConv = isapprox(μOld, μ, atol=μtol) |> (x -> sum(x) == length(x))
    EConv = isapprox(EOld, E, atol=Etol) |> (x -> sum(x) == length(x))

  end

end