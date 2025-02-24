"""
SCME/f potential 

This calls a ASE calculator of the SCME/f potential (not great)
I still want to make a better wrapper.

SCME Paper:
https://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp52097h

SCME/f Paper:
https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00598
"""

struct _SCMEf_PotVars
  lat
end

function SCMEf(cell::JMD.MyCell)
  lat = reshape(cell.lattice, 9)

  _SCMEf_PotVars(lat[[1,5,9]])
end

function SCMFf(F, G, y0, p)

  E = 0.0

  if G != nothing
    e, f = py"scmef_get_energy_and_forces"(y0, p.potVars.lat)
    E    = e
    G   .= f 
  else
    E    = py"scmef_get_energy"(y0, p.potVars.lat)
  end

  if F != nothing
    return E
  end

end

function SCMEf(F, G, cell::JMD.MyCell, lat)
end