
function pbcVec!(
  rbuf::AbstractVector, ri::AV, rj::AV, lat::AbstractMatrix
) where {AV<:AbstractVector}
  @. rbuf = rj - ri
  d       = norm(rbuf)

  if iszero(lat)
    return d
  end

  v = zeros(3)
  r = zeros(3)

  for i = -1:1
    for j = -1:1
      for k = -1:1
        v .= (lat[1, :] * i) + (lat[2, :] * j) + (lat[3, :] * k)
        r .= (rj .+ v) .- ri
        di = norm(r)

        if di < d
          d     = di
          rbuf .= r
        end
      end
    end
  end

  d
end
