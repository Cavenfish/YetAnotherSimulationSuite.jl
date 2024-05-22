"""
CO-CO Potential from van Hemert et al. 2015
"""

struct _MvHffCO_PotVars{F<:Float64} <: PotVars
  D::F
  a::F
  req::F
  Acc::F
  Aco::F
  Aoo::F
  Bcc::F
  Bco::F
  Boo::F
  Ccc::F
  Cco::F
  Coo::F
  Qc::F
  Qo::F
  αc::F
  αo::F
end

function MvHffCO(bdys::Vector{Atom}) = _MvHffCO_PotVars (
  11.230139012256362,
  2.3281,
  1.1282058129221093,
  361.367206403597,
  1516.76265699823,
  6370.185468304371,
  2.8345891887553925,
  3.5432364859442407,
  4.2518837831330885,
  33.44955570065988,
  15.189133724736946,
  10.546349734130885,
  3.843702939952312,
  2.131611069944055
)

function MvHffCO(F, G, y0, p)

  # initialize things
  P      = p.potVars
  E      = 0.0
  u      = Vector[]
  forces = Vector[]
  for i in 1:3:length(y0)
    push!(u, y0[i:i+2])
    push!(forces, [0.0, 0.0, 0.0])
  end

  for mol in p.mols
    E += _Morse!(forces, u, mol[1], mol[2], P.D, P.a, P.req)
  end

  for par in p.pars
    c1, o1 = p.pars[1]
    c2, o2 = p.pars[2]

    #C--C
    E += _shortDisp!(forces, u, c1, c2, P.Acc, P.Bcc)
    E += _longDisp!( forces, u, c1, c2, P.Ccc)

    #C--O
    E += _shortDisp!(forces, u, c1, o2, P.Aco, P.Bco)
    E += _longDisp!( forces, u, c1, o2, P.Cco)

    #O--C
    E += _shortDisp!(forces, u, o1, c2, P.Aco, P.Bco)
    E += _longDisp!( forces, u, o1, c2, P.Cco)

    #O--O
    E += _shortDisp!(forces, u, o1, o2, P.Aoo, P.Boo)
    E += _longDisp!( forces, u, o1, o2, P.Coo)

    #Special Electrostatics
    E += _electroMvH!(forces, u, par, P.Qc, P.Qo, P.αc, P.αo)

  end

  if G != nothing
    tmp = [j for i in forces for j in i]
    for i in 1:length(G)
      G[i] = -tmp[i]
    end
  end

  if F != nothing
    return E
  end

end

function _electroMvH!(F, u, par, Qc, Qo, αc, αo, req)

  #Atom indicies
  C1, O1 = pars[1]
  C2, O2 = pars[2]

  #Atom locations
  c1, o1 = u[pars[1]]
  c2, o2 = u[pars[2]]

  #Bond lengths
  R1 = c1 - o1
  R2 = c2 - o2
  r1 = norm(R1)
  r2 = norm(R2)

  #CM weights
  wo = 0.57135
  wc = 0.42865

  #X locations
  x1 = (wc*c1) + (wo*o1)
  x2 = (wc*c2) + (wo*o2)

  #Variable charges
  Qc1 = Qc * exp(-αc*(r1 - req))
  Qo1 = Qo * exp(-αo*(r1 - req))
  Qx1 = -(Qc1 + Qo1)
  Qc2 = Qc * exp(-αc*(r2 - req))
  Qo2 = Qo * exp(-αo*(r2 - req))
  Qx2 = -(Qc2 + Qo2)

  #The additional contribution to the forces resulting from
  #the position dependent charges is: -α * Q * rvec / r

  #C--C
  ϵ, F   = _Coulomb(c1, c2, Qc1, Qc2)
  E     += ϵ
  F_Q1   = αC * ϵ * R1 / r1
  F_Q2   = αC * ϵ * R2 / r2
  F[C1] -= (F + F_Q1)
  F[O1] += F_Q1
  F[C2] +=  F - F_Q2
  F[O2] += F_Q2

# C-X
  ϵ, F   = V_coul(rX2i0..., QC1*QX2)
  E     += ϵ
  F_Q1   = αC * ϵ * R1 / r1
  F_Q2   = - (αC * QC2 + αO * QO2) * ϵ/QX2 * R2 / r2
  F[C1] -= (F + F_Q1)
  F[O1] += F_Q1
  F[C2] += (wC * F) - F_Q2
  F[O2] += (wO * F) + F_Q2

  # C-O
  ϵ, F   = V_coul(rj1i0..., QC1*QO2)
  E     += ϵ
  F_Q1   = αC * ϵ * R1 / r1
  F_Q2   = αO * ϵ * R2 / r2
  F[C1] -= (F + F_Q1)
  F[O1] += F_Q1
  F[C2] -= F_Q2
  F[O2] += (F + F_Q2)

  # X-C
  ϵ, F   = V_coul(rj0X1..., QX1*QC2)
  E     += ϵ
  F_Q1   = - (αC * QC1 + αO * QO1) * ϵ/QX1 * R1 / r1
  F_Q2   = αC * ϵ * R2 / r2
  F[C1] -= ((wC * F) + F_Q1)
  F[O1] += F_Q1 - (wO * F)
  F[C2] += (F - F_Q2)
  F[O2] += F_Q2

  # X-X
  ϵ, F   = V_coul(rX2X1..., QX1*QX2)
  E     += ϵ
  F_Q1   = - (αC * QC1 + αO * QO1) * ϵ/QX1 * R1 / r1
  F_Q2   = - (αC * QC2 + αO * QO2) * ϵ/QX2 * R2 / r2
  F[C1] += ((-wC * F) - F_Q1)
  F[O1] += ((-wO * F) + F_Q1)
  F[C2] += (( wC * F) - F_Q2)
  F[O2] += (( wO * F) + F_Q2)


  # X-O
  ϵ, F   = V_coul(rj1X1..., QX1*QO2)
  E     += ϵ
  F_Q1   = - (αC * QC1 + αO * QO1) * ϵ/QX1 * R1 / r1
  F_Q2   = αO * ϵ * R2 / r2
  F[C1] -= ((wC * F) + F_Q1)
  F[O1] += (F_Q1 - (wO * F))
  F[C2] -= F_Q2
  F[O2] += (F + F_Q2)

  # O-C
  ϵ, F   = V_coul(rj0i1..., QO1*QC2)
  E     += ϵ
  F_Q1   = αO * E * R1 / r1
  F_Q2   = αC * E * R2 / r2
  F[C1] -= F_Q1
  F[O1] += (F_Q1 - F)
  F[C2] += (F - F_Q2)
  F[O2] += F_Q2

  # O-X
  ϵ, F   = V_coul(rX2i1..., QO1*QX2)
  E     += ϵ
  F_Q1   = αO * E * R1 / r1
  F_Q2   = - (αC * QC2 + αO * QO2) * E/QX2 * R2 / r2
  F[C1] -= F_Q1
  F[O1] += (F_Q1 - F)
  F[C2] += ((wC * F) - F_Q2)
  F[O2] += ((wO * F) + F_Q2)

  # O-O
  ϵ, F   = V_coul(rj1i1..., QO1*QO2)
  E     += ϵ
  F_Q1   = αO * E * R1 / r1
  F_Q2   = αO * E * R2 / r2
  F[C1] -= F_Q1
  F[O1] += (F_Q1 - F)
  F[C2] -= F_Q2
  F[O2] += (F + F_Q2)
  
  E
end