
function getVibCoup(EoM, bdys, ϵ, vec; r0=0.0)

  function bar(r)
    #make copy to prevent repeated changes
    tmp = deepcopy(bdys)

    #vec length (l) = 3N ; always 2 mols => l/6 = n
    #where n is number of atoms per mol
    n = div(length(vec), 6)

    #move along vec
    for i in 1:3:length(vec)
      j::UInt32 = (i+2)/3
      tmp[j].r += r * vec[i:i+2]
    end

    #get mol1 and mol2
    mol1 = tmp[1:n]
    mol2 = tmp[n+1:end]

    molsE = getPotEnergy(EoM, mol1) + getPotEnergy(EoM, mol2)
    E     = getPotEnergy(EoM, tmp) - molsE

    log(E + 2ϵ)
  end

  β = central_fdm(5, 1)(bar, r0)

  β
end