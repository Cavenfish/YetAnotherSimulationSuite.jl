"""
    diffDotSqrt(v2, v1)

Compute the vector difference and its Euclidean norm between two vectors.

# Arguments
- `v2`, `v1`: Input vectors.

# Returns
- Tuple: (Euclidean distance, difference vector).
"""
function diffDotSqrt(v2, v1)
  rvec = v2 - v1
  r    = sqrt(dot(rvec, rvec))
  (r, rvec)
end

"""
    getAngle(r1, r2)

Compute the angle (in radians) between two vectors.

# Arguments
- `r1`, `r2`: Input vectors.

# Returns
- Angle in radians.
"""
function getAngle(r1, r2)
  x  = (dot(r1, r2)) / (norm(r1) * norm(r2))
  x  = round(x, digits=15)
  
  acos(x)
end

