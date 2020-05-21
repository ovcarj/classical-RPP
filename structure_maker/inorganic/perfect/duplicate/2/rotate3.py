import numpy as np
from ase.io import read

atoms = read('2n3_final.traj')

cellx = atoms.get_cell()[0][0]
celly = atoms.get_cell()[1][1]
cellz = atoms.get_cell()[2][2]

up = atoms.copy()

up.positions[:,2] += 29.
#up.euler_rotate(phi=90., center='COM')
up.positions[:,0] += 0.5 * cellx
up.positions[:,1] += 0.5 * celly

atomss = atoms + up
atomss.set_cell(atoms.get_cell())

atomss.cell[2][2] += 28.

atomss.write('2n3_offset.traj')
