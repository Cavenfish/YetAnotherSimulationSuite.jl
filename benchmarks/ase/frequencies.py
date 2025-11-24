from ase.calculators.lj import LennardJones
from ase.io import read
from ase.vibrations import Vibrations
import sys

i = sys.argv[1]

bdys = read(f"../xyzfiles/{i}au_opt.xyz")

calc = LennardJones(sigma=2.95, epsilon=0.2297, ro=19.0, rc=20.0)

bdys.calc = calc

vib = Vibrations(bdys)
vib.run()

