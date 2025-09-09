
function pbcVec!(
  rbuf::AV3D, ri::AV3D, rj::AV3D, lat::AbstractMatrix
) where AV3D
  @. rbuf = rj - ri
  d       = norm(rbuf)

  if iszero(lat)
    return d
  end

  v  = zeros(3)
  rp = zeros(3)
  rn = zeros(3)

  for k = 1:3
    v  .= lat[k,:]
    rp .= (u[j] .+ v) .- u[i]
    rn .= (u[j] .- v) .- u[i]
    dp  = norm(rp)
    dn  = norm(rn)

    if dp < d && dp < dn
      d     = dp
      rbuf .= rp
    elseif dn < d
      d     = dn
      rbuf .= rn
    end
  end

  d
end
