
struct _PBC_Cache{AV3D<:AbstractVector}
  v::AV3D
  r::AV3D
end

const _pbc_buf = _PBC_Cache(MVector{3}(zeros(3)), MVector{3}(zeros(3)))

function pbcVec!(
  rbuf::AbstractVector, ri::AV, rj::AV, lat::AbstractMatrix;
  buf=_pbc_buf
) where {AV<:AbstractVector}
  @. rbuf = rj - ri
  d       = norm(rbuf)

  if iszero(lat)
    return d
  end

  for i = -1:1
    for j = -1:1
      for k = -1:1
        buf.v .= (lat[1, :] * i) + (lat[2, :] * j) + (lat[3, :] * k)
        buf.r .= (rj .+ buf.v) .- ri
        di     = norm(buf.r)

        if di < d
          d     = di
          rbuf .= buf.r
        end
      end
    end
  end

  d
end
