using StaticArrays

function calcMorse(r, diff)
  D      = 11.230139012256362
  a      = 2.626624 # 2.3281 * 1.1282
  r0     = 1.1282058129221093
  preF   = 2 * D * a / r0
  expf   = exp(a * (1 - r / r0))
  energy = D * expf * (expf - 2)
  F      = preF * expf * (expf - 1) * diff / r

  return E,F
end

function V_disp(diff, Cij)
  r2 = diff'diff
  E  = Cij / r2^3
  F  = 6.0 * E * diff / r2
  return E,F
end

function calcDisp(rj0i0, rj1i1, rj0i1, rj1i0)
  Ccc = -33.44955570065988
  Coo = -10.546349734130885
  Cco = -15.189133724736946
  Coc = -15.189133724736946

  # C-C
  E, Fcc  = V_disp(rj0i0..., Ccc)
  energy  = E

  # O-O
  E, Foo  = V_disp(rj1i1..., Coo)
  energy += E

  # C-O
  E, Fco  = V_disp(rj0i1..., Cco)
  energy += E

  # O-C
  E, Foc  = V_disp(rj1i0..., Coc)
  energy += E

  return energy, Fcc, Foo, Fco, Foc
end

function V_exch(diff, Aij, Bij)
  r = √(diff'diff)
  E = Aij * exp(-Bij * r)
  F = Bij * E * diff / r
  return E,F
end

function calcExch(rj0i0, rj1i1, rj0i1, rj1i0)
  Acc = 361.367206403597
  Aoo = 6370.185468304371
  Aco = 1516.76265699823
  Aoc = 1516.76265699823
  Bcc = 2.8345891887553925
  Boo = 4.2518837831330885
  Bco = 3.5432364859442407
  Boc = 3.5432364859442407

  # C-C
  E, Fcc  = V_exch(rj0i0..., Acc, Bcc)
  energy  = E

  # O-O
  E, Foo  = V_exch(rj1i1..., Aoo, Boo)
  energy += E

  # C-O
  E, Fco  = V_exch(rj0i1..., Aco, Bco)
  energy += E

  # O-C
  E, Foc  = V_exch(rj1i0..., Aoc, Boc)
  energy += E

  return energy, Fcc, Foo, Fco, Foc
end

function V_coul(r, diff, Qij)
  E = Qij / r
  F = Qij * diff / r**3 	# = Qij / r**2 * diff / r
  return E, F
end

function calcCoul(rC1, rO1, rC2, rO2, rj0i0, rj1i1, rj0i1, rj1i0, ri1i0)
  r0          =  1.1282058129221093
  Qc          = -1.7835026375774934
  Qo          = -2.333732174702465
  alphaC      =  3.843702939952312
  alphaO      =  2.131611069944055
  wO          =  0.57135
  wC          =  0.42865

  ### CO1 ###
# atom positions
  r1, rCO1 = ri1i0

  # gradient of r1 with respect to rC1: - rCO1 / r1
# gradient of r1 with respect to rO1: rCO1 / r1

# center of mass
  rX1 = (wC*rC1) + (wO*rO1)

# rCO1-dependent charges:
  QC1 = Qc * exp(-alphaC*(r1 - r0))
  QO1 = Qo * exp(-alphaO*(r1 - r0))
  QX1 = - (QC1 + QO1)

  ### CO2 ###
  # atom positions
  r2, rCO2 = diffDotSqrt(rO2, rC2)

  # gradient of r2 with respect to rC2: - rCO2 / r2
  # gradient of r2 with respect to rO2: rCO2 / r2

  # center of mass
  rX2 = (wC*rC2) + (wO*rO2)

  # rCO2-dependent charges:
  QC2 = Qc * exp(-alphaC*(r2 - r0))
  QO2 = Qo * exp(-alphaO*(r2 - r0))
  QX2 = - (QC2 + QO2)

  # for force contributions resulting from bondlength-dependent charges:
  # nabla_r Q(r) = -alpha * Q(r) * nabla_r r

  # rX i/j preps
  rX2i0 = diffDotSqrt(rX2,rC1)
  rX2i1 = diffDotSqrt(rX2,rO1)
  rj0X1 = diffDotSqrt(rC2,rX1)
  rj1X1 = diffDotSqrt(rO2,rX1)
  rX2X1 = diffDotSqrt(rX2,rX1)

# C-C
  E, F    = V_coul(*rj0i0, QC1*QC2)
  energy  = E
  F_Q1    = alphaC * E * rCO1 / r1
  F_Q2    = alphaC * E * rCO2 / r2
  Fi0     = - (F + F_Q1)
  Fj0     =    F - F_Q2
  Fi1     = F_Q1
  Fj1     = F_Q2

# C-X
  E, F    = V_coul(*rX2i0, QC1*QX2)
  energy += E
  F_Q1    = alphaC * E * rCO1 / r1
  F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
  Fi0    -= (F + F_Q1)
  Fj0    += ((wC * F) - F_Q2)
  Fi1    += F_Q1
  Fj1    += (wO * F) + F_Q2

  # C-O
  E, F    = V_coul(*rj1i0, QC1*QO2)
  energy += E
  F_Q1    = alphaC * E * rCO1 / r1
  F_Q2    = alphaO * E * rCO2 / r2
  Fi0    -= (F + F_Q1)
  Fj1    += (F + F_Q2)
  Fi1    += F_Q1
  Fj0    -= F_Q2

  # X-C
  E, F    = V_coul(*rj0X1, QX1*QC2)
  energy += E
  F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
  F_Q2    = alphaC * E * rCO2 / r2
  Fi0    -= ((wC * F) + F_Q1)
  Fi1    += F_Q1 - (wO * F)
  Fj0    += (F - F_Q2)
  Fj1    += F_Q2

  # X-X
  E, F    = V_coul(*rX2X1, QX1*QX2)
  energy += E
  F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
  F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
  Fi0    += ((-wC * F) - F_Q1)
  Fi1    += ((-wO * F) + F_Q1)
  Fj0    += (( wC * F) - F_Q2)
  Fj1    += (( wO * F) + F_Q2)


  # X-O
  E, F    = V_coul(*rj1X1, QX1*QO2)
  energy += E
  F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
  F_Q2    = alphaO * E * rCO2 / r2
  Fi0    -= ((wC * F) + F_Q1)
  Fi1    += (F_Q1 - (wO * F))
  Fj0    -= F_Q2
  Fj1    += (F + F_Q2)

  # O-C
  E, F    = V_coul(*rj0i1, QO1*QC2)
  energy += E
  F_Q1    = alphaO * E * rCO1 / r1
  F_Q2    = alphaC * E * rCO2 / r2
  Fi0    -= F_Q1
  Fi1    += (F_Q1 - F)
  Fj0    += (F - F_Q2)
  Fj1    += F_Q2

  # O-X
  E, F    = V_coul(*rX2i1, QO1*QX2)
  energy += E
  F_Q1    = alphaO * E * rCO1 / r1
  F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
  Fi0    -= F_Q1
  Fi1    += (F_Q1 - F)
  Fj0    += ((wC * F) - F_Q2)
  Fj1    += ((wO * F) + F_Q2)

  # O-O
  E, F    = V_coul(*rj1i1, QO1*QO2)
  energy += E
  F_Q1    = alphaO * E * rCO1 / r1
  F_Q2    = alphaO * E * rCO2 / r2
  Fi0    -= F_Q1
  Fi1    += (F_Q1 - F)
  Fj0    -= F_Q2
  Fj1    += (F + F_Q2)

  return energy, Fi0, Fi1, Fj0, Fj1
end

function diffDotSqrt(v2, v1)
  diff = v2 - v1
  r    = sqrt(dot(diff, diff))
  return (r, diff)
end

function COCOdyn(dv, v, u, p, t)
  epsilon = 11.230139012256362
  N       = len(positions)
  MorseF  = zeros_like(positions)
  ExchF   = zeros_like(positions)
  DispF   = zeros_like(positions)
  CoulF   = zeros_like(positions)
  MorseE  = ExchE = DispE = CoulE = 0.0

  for i in range(0, N//2, 1):
    posi0 = positions[2*i]
    posi1 = positions[2*i+1]
    ri1i0 = diffDotSqrt(posi1, posi0)

    #Calculate Morse
    E, F           = calculate_Morse(*ri1i0)
    MorseE        += E
    MorseF[2*i]   -= F
    MorseF[2*i+1] += F

    for j in range(i+1, N//2, 1):
      posj0 = positions[2*j]
      posj1 = positions[2*j+1]

      rj0i0 = diffDotSqrt(posj0, posi0)
      rj1i1 = diffDotSqrt(posj1, posi1)
      rj0i1 = diffDotSqrt(posj0, posi1)
      rj1i0 = diffDotSqrt(posj1, posi0)

      #Calculate Dispersion
      E, Fcc, Foo, Fco, Foc = calculate_Disp(rj0i0, rj1i1, rj0i1, rj1i0)
      DispE                += E
      DispF[2*i]           -= (Fcc + Foc)
      DispF[2*i+1]         -= (Foo + Fco)
      DispF[2*j]           += (Fcc + Fco)
      DispF[2*j+1]         += (Foo + Foc)

      #Calculate Exchange
      E, Fcc, Foo, Fco, Foc = calculate_Exch(rj0i0, rj1i1, rj0i1, rj1i0)
      ExchE                += E
      ExchF[2*i]           -= (Fcc + Foc)
      ExchF[2*i+1]         -= (Foo + Fco)
      ExchF[2*j]           += (Fcc + Fco)
      ExchF[2*j+1]         += (Foo + Foc)

      #Calculate Coulomb
      E, Fi0, Fi1, Fj0, Fj1 = calculate_Coul(posi0, posi1, posj0, posj1,
                                             rj0i0, rj1i1, rj0i1, rj1i0,
                                             ri1i0)
      CoulE                += E
      CoulF[2*i]           += Fi0
      CoulF[2*i+1]         += Fi1
      CoulF[2*j]           += Fj0
      CoulF[2*j+1]         += Fj1

  # In order to normalize the minimal morse potential to zero,
  # the following energy will be added to change the zero point
  # of all intramolecular interactions.
  MorseE += epsilon * N / 2.0

  Morse = (MorseE, MorseF)
  Exch  = ( ExchE,  ExchF)
  Disp  = ( DispE,  DispF)
  Coul  = ( CoulE,  CoulF)
  return Morse, Exch, Disp, Coul
end
