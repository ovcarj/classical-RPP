import numpy as np
from ase.io import read

atoms = read('2n2_final.traj')

N_indices = np.where(atoms.symbols == 'N')[0][0:8]
z0 = np.average(atoms.positions[N_indices][:,2])
above_N = N_indices[np.where(atoms.positions[:,2][N_indices] > z0)[0]]


cellx = atoms.get_cell()[0][0]
celly = atoms.get_cell()[1][1]
cellz = atoms.get_cell()[2][2]

up = atoms.copy()

up.positions[:,2] += 24.
#up.euler_rotate(phi=90., center='COM')
up.positions[:,0] += 0.25 * cellx
up.positions[:,1] += 0.25 * celly

atomss = atoms + up
atomss.set_cell(atoms.get_cell())

atomss.cell[2][2] += 24.

atomss.write('2n2_offset_super.traj')
