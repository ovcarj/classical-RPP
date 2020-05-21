import numpy as np
from ase.io import read
from ase import neighborlist
from scipy import sparse
from ase import Atoms
from ase.build import sort

frame = read('1n3_MA.traj')
MA_atoms = Atoms()

#up = atoms.copy()

#up.positions[:,2] += 22.

#atomss = atoms + up
#atomss.set_cell(atoms.get_cell())

#atomss.cell[2][2] += 22.

#atomss.write('1n2_up.traj')

symbols = np.asarray(frame.get_chemical_symbols(), dtype='str')

Pb_indices = np.flatnonzero(symbols == 'Pb')
Br_indices = np.flatnonzero(symbols == 'Br')
N_indices = np.flatnonzero(symbols == 'N')
C_indices = np.flatnonzero(symbols == 'C')
H_indices = np.flatnonzero(symbols == 'H')

organic_indices = np.concatenate((N_indices, C_indices, H_indices))

organic = frame[organic_indices]

inorganic_indices = np.concatenate((Pb_indices, Br_indices))

inorganic = frame[inorganic_indices]

cutOff = neighborlist.natural_cutoffs(organic)
neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
neighborList.update(organic)

matrix = neighborList.get_connectivity_matrix()

n_components, component_list = sparse.csgraph.connected_components(matrix)

MA=np.empty((0), dtype='int')
MA_counter = 0

for i in range(n_components):
    molIdx=i
    molIdxs = [ j for j in range(len(component_list)) if component_list[j] == molIdx ]
    if(len(molIdxs) == 8):
        MA = np.append(MA, molIdxs, axis=0)
        MA_counter += 1

if(len(MA) > 0):
    MA = np.split(MA, MA_counter)

for i in range(len(MA)):
    MA_atoms += sort(organic[MA[i]])

MA_atoms_up = MA_atoms.copy()
MA_atoms_up.positions[:,2] += 27.

inorganic_up = inorganic.copy()
inorganic_up.positions[:,2] += 27.

atoms = MA_atoms + MA_atoms_up + inorganic + inorganic_up

atoms.set_cell(frame.get_cell())
atoms.cell[2][2] += 27.

atoms.write('1n3_up.traj')
