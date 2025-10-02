"""
# Simple Point Charge (SPC) Models

## SPC/F

Based on:
Toukan, Kahled, and Aneesur Rahman. "Molecular-dynamics study of atomic
motions in water." Physical Review B 31.5 (1985): 2643.

Link:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.31.2643

## SPC/Fd

Based on:
Dang, Liem X., and B. Montgomery Pettitt. "Simple intramolecular model
potentials for water." Journal of physical chemistry 91.12 (1987): 3349-3354.

Link:
https://doi.org/10.1021/j100296a048

## SPC/Fw

Based on:
Wu, Yujie, Harald L. Tepper, and Gregory A. Voth. "Flexible simple
point-charge water model with improved liquid-state properties." The Journal
of chemical physics 124.2 (2006).

Link:
https://doi.org/10.1063/1.2136877
"""
SPC(name::String; constraints=nothing) = Calculator(
  x -> SPC(x, name), 
  u"eV", u"eV/angstrom", u"angstrom * u^0.5 * eV^-0.5";
  EF=SPC_ef!, constraints=constraints
)

struct _SPCF_PotVars{
  F<:Float64, AV3D<:AbstractVector, AV2D<:AbstractVector
} <: PotVars
  Kb::F
  req::F
  Kθ::F
  θeq::F
  σ::F
  ϵ::F
  Qo::F
  Qh::F
  rs::F
  rc::F
  S::AV2D
  Fbuf::AV3D
  rbuf::AV3D
end

function SPC(bdys::Union{Vector{MyAtoms}, MyCell}, name::String) 
  
  if name == "SPC/F"
    return _SPCF_PotVars(
      48.05913,
      1.0,
      3.97,
      1.910611,
      3.145, 
      0.007,
      -2.959855,
      0.5 * 2.959855,
      9.0,
      10.0,
      MVector{2}(zeros(2)),
      MVector{3}(zeros(3)),
      MVector{3}(zeros(3))
    )
  end

  if name == "SPC/Fd"
    return _SPCF_PotVars(
      45.71445, # g
      1.0, # g
      3.2913, # g
      1.911135, # g
      3.165492, # g
      0.0067398, # g
      -3.1116428, # g
      0.5 * 3.1116428, # g
      9.0,
      10.0,
      MVector{2}(zeros(2)),
      MVector{3}(zeros(3)),
      MVector{3}(zeros(3))
    )
  end

  if name == "SPC/Fw"
    return _SPCF_PotVars(
      45.9296231, # g
      1.012, # g
      3.2913, # g
      1.976410, # g
      3.165492, # g
      0.0067398, # g
      -3.1116428, # g
      0.5 * 3.1116428, # g
      9.0,
      10.0,
      MVector{2}(zeros(2)),
      MVector{3}(zeros(3)),
      MVector{3}(zeros(3))
    )
  end

  ArgumentError("$(name) is not an available SPC potential") |> throw

end

function SPC_ef!(F, u, p)
  E   = 0.0
  P   = p.potVars
  lat = isa(p, optVars) ? p.lattice : p.ensemble.lattice

  for i = 1:length(p.mols)
    o1, h1, h2 = p.mols[i]

    E += harmonicBond!(F, u, lat, o1, h1, P.Kb, P.req, P.rc)
    E += harmonicBond!(F, u, lat, o1, h2, P.Kb, P.req, P.rc)
    E += harmonicBondAngle!(F, u, lat, h1, o1, h2, P.Kθ, P.θeq, P.rc)

    for j = i+1:length(p.mols)
      o2, h3, h4 = p.mols[j]

      doo = pbcVec!(P.rbuf, u[o1], u[o2], P.rc, lat)
      switchLR!(P.S, doo, P.rs, P.rc)

      if iszero(P.S)
        continue
      end

      e  = 0.0

      # O-O Interactions
      e += vdw!(F, u, lat, o1, o2, P.ϵ, P.σ, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, o1, o2, P.Qo, P.Qo, P.rc; S=P.S[1])

      # H-H Coulomb
      e += coulomb!(F, u, lat, h1, h3, P.Qh, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, h1, h4, P.Qh, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, h2, h3, P.Qh, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, h2, h4, P.Qh, P.Qh, P.rc; S=P.S[1])

      # H-O Coulomb
      e += coulomb!(F, u, lat, o1, h3, P.Qo, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, o1, h4, P.Qo, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, o2, h1, P.Qo, P.Qh, P.rc; S=P.S[1])
      e += coulomb!(F, u, lat, o2, h2, P.Qo, P.Qh, P.rc; S=P.S[1])

      E += P.S[1] * e

      if P.S[2] != 0.0
        P.rbuf ./= doo
        P.Fbuf  .= -P.S[2] * e .* P.rbuf
        F[o1]  .-= P.Fbuf
        F[o2]  .+= P.Fbuf
      end
    end
  end

  E
end