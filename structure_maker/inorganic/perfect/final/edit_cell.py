import numpy as np
import sys
from ase.io import read

atoms = read(str(sys.argv[1]))

atoms.cell[2][2] = np.max(atoms.positions[:,2])

atoms.write(str(sys.argv[1]))
