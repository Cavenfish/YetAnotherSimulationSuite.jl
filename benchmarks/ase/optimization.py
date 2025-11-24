from ase.calculators.lj import LennardJones
from ase.io import read
from ase.optimize import LBFGS  
import sys

i = sys.argv[1]
j = int(sys.argv[2])

bdys = read(f"../xyzfiles/{i}au.xyz")

calc = LennardJones(sigma=2.95, epsilon=0.2297, ro=19.0, rc=20.0)

bdys.calc = calc

dyn = LBFGS(bdys, logfile=None)
dyn.run(fmax=1e-18, steps=j) 