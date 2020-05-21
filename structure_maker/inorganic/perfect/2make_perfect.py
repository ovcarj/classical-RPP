import numpy as np

from ase.io import read, write
from ase import Atoms

#dft = read('n1.traj')

Pb_Br = 3.040 #Pb-Br z-distance

cell_x = 3. * Pb_Br
cell_y = 3. * Pb_Br
cell_z = 2. * Pb_Br

no_Pb = str(4)
no_Br = str(16)

atoms = Atoms('Br' + no_Br + 'Pb' + no_Pb, [
 (0. * Pb_Br, 3. * Pb_Br, 3. * Pb_Br),
 (0. * Pb_Br, 1. * Pb_Br, 3. * Pb_Br),
 (1. * Pb_Br, 2. * Pb_Br, 3. * Pb_Br),
 (1. * Pb_Br, 0. * Pb_Br, 3. * Pb_Br),
 (2. * Pb_Br, 3. * Pb_Br, 3. * Pb_Br),
 (2. * Pb_Br, 1. * Pb_Br, 3. * Pb_Br),
 (3. * Pb_Br, 2. * Pb_Br, 3. * Pb_Br),
 (3. * Pb_Br, 0. * Pb_Br, 3. * Pb_Br),
 (1. * Pb_Br, 3. * Pb_Br, 4. * Pb_Br),
 (1. * Pb_Br, 1. * Pb_Br, 4. * Pb_Br),
 (3. * Pb_Br, 3. * Pb_Br, 4. * Pb_Br),
 (3. * Pb_Br, 1. * Pb_Br, 4. * Pb_Br),
 (1. * Pb_Br, 3. * Pb_Br, 2. * Pb_Br),
 (1. * Pb_Br, 1. * Pb_Br, 2. * Pb_Br),
 (3. * Pb_Br, 3. * Pb_Br, 2. * Pb_Br),
 (3. * Pb_Br, 1. * Pb_Br, 2. * Pb_Br),
 (1. * Pb_Br, 3. * Pb_Br, 3. * Pb_Br),
 (1. * Pb_Br, 1. * Pb_Br, 3. * Pb_Br),
 (3. * Pb_Br, 3. * Pb_Br, 3. * Pb_Br),
 (3. * Pb_Br, 1. * Pb_Br, 3. * Pb_Br)])



atoms.set_cell([cell_x, cell_y, cell_z])

atoms.write('2n1.traj')
