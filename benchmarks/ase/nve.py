from ase.calculators.lj import LennardJones
from ase.io import read
import sys
from ase.md.verlet import VelocityVerlet
from ase import units

i = sys.argv[1]
j = int(sys.argv[2]) * 1000

bdys = read(f"../xyzfiles/{i}au.xyz")

calc = LennardJones(sigma=2.95, epsilon=0.2297, ro=19.0, rc=20.0)

bdys.calc = calc

dyn = VelocityVerlet(bdys, 1 * units.fs)
dyn.run(j)