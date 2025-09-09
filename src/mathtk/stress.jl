"""
    getNumericalStress(calc::MyCalc, cell::MyCell; eps=1e-6)

Compute the numerical stress tensor for a cell using finite differences.

# Arguments
- `calc`: Calculator object (`MyCalc`).
- `cell`: Cell object (`MyCell`).
- `eps`: Strain increment for finite difference (default: 1e-6).

# Returns
- 3x3 stress tensor matrix.
"""
function getNumericalStress(calc::MyCalc, cell::MyCell; eps=1e-6)
  stress = zeros(3,3)
  vol    = getVolume(cell)
  tmp    = deepcopy(cell)

  for i = 1:3

    x = Matrix(1.0I, 3,3)

    x[i,i]       = 1.0 + eps
    tmp.lattice .= cell.lattice * x
    eplus        = getPotEnergy(calc, tmp)

    x[i,i]       = 1.0 - eps
    tmp.lattice .= cell.lattice * x
    eminus       = getPotEnergy(calc, tmp)

    stress[i,i]  = (eplus - eminus) / (2 * eps * vol)
    x[i,i]       = 1.0

    i+1 > 3 ? j =  (i+1)%3 : j = i+1

    x[i,j] = x[j,i] = +0.5 * eps
    tmp.lattice    .= cell.lattice * x
    eplus           = getPotEnergy(calc, tmp)

    x[i,j] = x[j,i] = -0.5 * eps
    tmp.lattice    .= cell.lattice * x
    eminus          = getPotEnergy(calc, tmp)

    stress[i,j] = stress[j,i] = (eplus - eminus) / (2 * eps * vol)
  end

  stress
end

"""
    getNumericalStressOrthogonal(calc::MyCalc, cell::MyCell; eps=1e-6)

Compute the numerical stress tensor for a cell using finite differences, only for diagonal terms (orthogonal cell).

# Arguments
- `calc`: Calculator object (`MyCalc`).
- `cell`: Cell object (`MyCell`).
- `eps`: Strain increment for finite difference (default: 1e-6).

# Returns
- 3x3 stress tensor matrix (only diagonal terms are nonzero).
"""
function getNumericalStressOrthogonal(calc::MyCalc, cell::MyCell; eps=1e-6)
  stress = zeros(3,3)
  vol    = getVolume(cell)
  tmp    = deepcopy(cell)

  for i = 1:3

    x = Matrix(1.0I, 3,3)

    x[i,i]       = 1.0 + eps
    tmp.lattice .= cell.lattice * x
    eplus        = getPotEnergy(calc, tmp)

    x[i,i]       = 1.0 - eps
    tmp.lattice .= cell.lattice * x
    eminus       = getPotEnergy(calc, tmp)

    stress[i,i]  = (eplus - eminus) / (2 * eps * vol)
    x[i,i]       = 1.0

  end

  stress
end