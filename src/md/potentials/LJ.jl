"""
Lennard Jones potential
"""
LJ(p::Dict; constraints=nothing) = Calculator(
  x -> build_lj_struct(x, p),
  u"eV", u"eV/angstrom", u"angstrom * u^0.5 * eV^-0.5";
  E=LJ_e, EF=LJ_ef!, constraints=constraints
)

struct _LJ_PotVars{
  F<: Float64, AV3D<:AbstractVector, AV2D<:AbstractVector
} <:PotVars
  ϵ::F
  σ::F
  rs::F
  rc::F
  S::AV2D
  Fbuf::AV3D
  rbuf::AV3D
end

function build_lj_struct(bdys::Union{MyCell, Vector{MyAtoms}}, p::Dict)
  _LJ_PotVars(
    p["epsilon"],
    p["sigma"],
    p["rs"],
    p["rc"],
    MVector{2}(zeros(2)),
    MVector{3}(zeros(3)),
    MVector{3}(zeros(3)),
  )
end

function LJ_e(u, p)
  E   = 0.0
  P   = p.potVars
  lat = isa(p, optVars) ? p.lattice : p.ensemble.lattice

  for i = 1:length(u)
    for j = i+1:length(u)

      d = pbcVec!(P.rbuf, u[i], u[j], P.rc, lat)
      switchLR!(P.S, d, P.rs, P.rc)

      if iszero(P.S)
        continue
      end

      a  = P.σ / d
      e  = 4*P.ϵ * ((a)^12 - (a)^6)
      E += P.S[1] * e
    end
  end
  
  E
end

function LJ_ef!(F, u, p)
  E   = 0.0
  P   = p.potVars
  lat = isa(p, optVars) ? p.lattice : p.ensemble.lattice

  for i = 1:length(u)
    for j = i+1:length(u)

      d = pbcVec!(P.rbuf, u[i], u[j], P.rc, lat)
      switchLR!(P.S, d, P.rs, P.rc)

      if iszero(P.S)
        continue
      end

      e = vdw!(F, u, lat, i, j, P.ϵ, P.σ, P.rc; S=P.S[1])

      E += P.S[1] * e

      if P.S[2] != 0.0
        P.rbuf ./= d
        P.Fbuf  .= -P.S[2] * e .* P.rbuf
        F[i]  .-= P.Fbuf
        F[j]  .+= P.Fbuf
      end
    end
  end
  
  E
end