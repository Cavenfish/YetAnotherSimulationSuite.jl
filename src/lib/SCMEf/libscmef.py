"""
SCME/f api for Julia 

load with PyCall using @pyinclude(FILE_NAME)

PyCall wrapping will make calls to SCME/f slower, but
making this wrapper is quicker than a proper c wrapper.
"""

def scmef_init(pos, cell, pbc=True, **kwargs):
    # SCMEf imports
    import pyscme
    from pyscme.parameters import parameter_H2O
    from pyscme.scme_calculator import SCMECalculator

    from ase import Atoms
    from ase.units import Bohr, Hartree

    # SCME params
    para_dict = {
        "te": 1.1045 / Bohr,
        "td": 7.5548 * Bohr,
        "Ar": 8149.63 / Hartree,
        "Br": -0.5515,
        "Cr": -3.4695 * Bohr,
        "r_Br": 1.0 / Bohr,
        "rc_Disp": 8.0 / Bohr,
        "rc_Core": 7.5 / Bohr,
        "rc_Elec": 9.0 / Bohr,
        "w_rc_Elec": 2.0 / Bohr,
        "w_rc_Core": 2.0 / Bohr,
        "w_rc_Disp": 2.0 / Bohr,
        "C6": 46.4430e0,
        "C8": 1141.7000e0,
        "C10": 33441.0000e0,
        "scf_convcrit": 1e-12,
        "scf_policy": pyscme.SCFPolicy.strict,
        "NC": [1,1,1],
        "dms": True,
        "qms": True
    }

    for (key, value) in kwargs.items():
        para_dict[key] = value

    # Define system
    symbols  = ["O", "H", "H"] * int(len(pos)/3)
    bdys     = Atoms(symbols, pos, cell=cell)
    bdys.pbc = pbc

    # Set SCME calc
    bdys.calc =  SCMECalculator(atoms=bdys, **para_dict)
    parameter_H2O.Assign_parameters_H20(bdys.calc.scme)

    return bdys

def scmef_get_energy(pos, cell, pbc=True, **kwargs):

    bdys = scmef_init(pos, cell, pbc=pbc, **kwargs)

    return bdys.get_potential_energy()

def scmef_get_energy_and_forces(pos, cell, pbc=True, **kwargs):

    bdys = scmef_init(pos, cell, pbc=pbc)

    E = bdys.get_potential_energy()
    F = bdys.get_forces()

    return E,F

def scmef_get_dipole(pos, cell, pbc=True, **kwargs):
    
    bdys = scmef_init(pos, cell, pbc=pbc)

    return bdys.get_dipole_moment()
    
def scmef_get_induced_dipoles(pos, cell, pbc=True, **kwargs):

    bdys = scmef_init(pos, cell, pbc=pbc)
    bdys.get_potential_energy()

    return bdys.calc.results["induced_dipoles"]

def scmef_get_constituent_energies(pos, cell, pbc=True, **kwargs):
    
    bdys = scmef_init(pos, cell, pbc=pbc)
    bdys.get_potential_energy()

    E_elec  = bdys.calc.results["energy_electrostatic"]
    E_disp  = bdys.calc.results["energy_dispersion"]
    E_core  = bdys.calc.results["energy_core"]
    E_intra = bdys.calc.results["energy_monomer"]

    return E_elec, E_disp, E_core, E_intra
