
function spreadForce!(
  F::AbstractVector, f::AbstractVector, 
  a::Vector{Int64}, wa::Vector{Float64}
) 

  for i = 1:length(a)
    F[a[i]] .+= f * wa[i]
  end

end