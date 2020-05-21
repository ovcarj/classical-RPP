import numpy as np
from ase.io import read

atoms = read('1n1.traj')

up = atoms.copy()

up.positions[:,2] += 16.
up.euler_rotate(phi=90., center='COM')

atomss = atoms + up
atomss.set_cell(atoms.get_cell())

atomss.cell[2][2] += 16.

atomss.write('1n1_rot.traj')
