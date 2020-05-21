import numpy as np

from ase.io import read, write
from ase import Atoms

#dft = read('n1.traj')

Pb_Br = 3.040 #Pb-Br z-distance

cell_x = 2. * np.sqrt(2) * Pb_Br
cell_y = 2. * np.sqrt(2) * Pb_Br
cell_z = 2. * Pb_Br

no_Pb = str(2)
no_Br = str(8)

atoms = Atoms('Br' + no_Br + 'Pb' + no_Pb, [
 (1./np.sqrt(2) * Pb_Br, 1./np.sqrt(2) * Pb_Br, 1.0 * Pb_Br),
 (3. * 1./np.sqrt(2) * Pb_Br, 3 * 1./np.sqrt(2) * Pb_Br, 1.0 * Pb_Br),
 (1./np.sqrt(2) * Pb_Br, 3. * 1./np.sqrt(2) * Pb_Br, 1.0 * Pb_Br),
 (3. * 1./np.sqrt(2) * Pb_Br, 1./np.sqrt(2) * Pb_Br, 1.0 * Pb_Br),
 (0.0 * Pb_Br, 0.0 * Pb_Br, 2. * Pb_Br),
 (2. * 1./np.sqrt(2) * Pb_Br, 2 * 1./np.sqrt(2) * Pb_Br, 2. * Pb_Br),
 (0.0 * Pb_Br, 0.0 * Pb_Br, 0.0 * Pb_Br),
 (2. * 1./np.sqrt(2) * Pb_Br, 2 * 1./np.sqrt(2) * Pb_Br, 0.0 * Pb_Br),
 (0.0 * Pb_Br, 0.0 * Pb_Br, 1.0 * Pb_Br),
 (2. * 1./np.sqrt(2) * Pb_Br, 2 * 1./np.sqrt(2) * Pb_Br, 1.0 * Pb_Br)])

atoms.set_cell([cell_x, cell_y, cell_z])

atoms.write('1n1.traj')
