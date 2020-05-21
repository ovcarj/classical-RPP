import numpy as np
from ase.io import read
from ase import neighborlist
from scipy import sparse
from ase import Atoms
from ase.build import sort

frame = read('2n3_final.traj')
MA_atoms = Atoms()
long_atoms = Atoms()

cell_x = frame.get_cell()[0][0]
cell_y = frame.get_cell()[1][1]

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

Long=np.empty((0), dtype='int')
Long_counter = 0

for i in range(n_components):
    molIdx=i
    molIdxs = [ j for j in range(len(component_list)) if component_list[j] == molIdx ]
    if(len(molIdxs) == 8):
        MA = np.append(MA, molIdxs, axis=0)
        MA_counter += 1
    if(len(molIdxs) > 8):
        Long = np.append(Long, molIdxs, axis=0)
        Long_counter += 1

if(len(MA) > 0):
    MA = np.split(MA, MA_counter)
    for i in range(len(MA)):
        MA_atoms += sort(organic[MA[i]])

Long = np.split(Long, Long_counter)

for i in range(len(Long)):
    long_atoms += sort(organic[Long[i]])

if(len(MA) > 0):
    MA_atoms_up = MA_atoms.copy()
    MA_atoms_up.positions[:,2] += 30.
    MA_atoms_up.positions[:,0] += 0.25 * cell_x
    MA_atoms_up.positions[:,1] += 0.25 * cell_y

#long_atoms_up = long_atoms.copy()
#long_atoms_up.positions[:,2] += 24.

z0 = np.average(long_atoms.positions[:,2])
long_above = long_atoms[np.where(long_atoms.positions[:,2] > z0)[0]]
long_below = long_atoms[np.where(long_atoms.positions[:,2] < z0)[0]]

long_above_up = long_above.copy()
long_below_up = long_below.copy()

long_above.positions[:,0] += 0.25 * cell_x
long_above.positions[:,1] += 0.25 * cell_y

#long_below.positions[:,0] += 0.25 * cell_x
#long_below.positions[:,1] += 0.25 * cell_y

long_above_up.positions[:,2] += 30.
long_below_up.positions[:,2] += 30.

long_below_up.positions[:,0] += 0.25 * cell_x
long_below_up.positions[:,1] += 0.25 * cell_y

inorganic_up = inorganic.copy()
inorganic_up.positions[:,2] += 30.
inorganic_up.positions[:,0] += 0.25 * cell_x
inorganic_up.positions[:,1] += 0.25 * cell_y


#atoms = MA_atoms + MA_atoms_up + inorganic + inorganic_up
atoms = long_below + long_above + long_below_up + long_above_up + MA_atoms + MA_atoms_up + inorganic + inorganic_up

atoms.set_cell(frame.get_cell())
#atoms.cell[0][0] += 2 * 3.040
#atoms.cell[1][1] += 2 * 3.040

atoms.cell[2][2] += 29.

atoms.write('2n3_offset_super.traj')
