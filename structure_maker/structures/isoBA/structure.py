#usr/bin/python

import numpy as np
import sys

from ase.io import read, write
from ase.io.trajectory import Trajectory, TrajectoryWriter
from ase import neighborlist
from ase import Atoms
from ase.build import sort

from scipy import sparse

def get_molecule_length(mol):
    N_index = np.where(mol.symbols == 'N')
    C_indices = np.where(mol.symbols == 'C')
    
    length = np.amax(mol.get_distances(N_index, C_indices))
    
    return length

def normalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

def get_rotation_matrix(a,b):
    scalar = np.dot(a,-b)
    cross = np.cross(a,-b)
    I = np.identity(3)
    V = np.array(([0., -cross[2], cross[1]],[cross[2], 0., -cross[0]],[-cross[1], cross[0], 0.]))
    VV = np.dot(V,V)
    
    R = I + V + VV * 1. / (1. + scalar)
    
    return(R)

ordered = Atoms()
upper_MA = Atoms()
bottom_MA = Atoms()

cell_x_factor = 0.0
cell_y_factor = 0.0

print('Enter unit cell type (1 or 2):')
print('X and Y cell lengths:')
print('1 = sqrt(2)')
print('2 = 2')

cell_type = str(input())

print('Enter number of inorganic layers (n):')

n = int(input())

print('Do you want the cell to be regular or offset? (enter r/o)')

ro = str(input())

if (ro == 'r'):
    reof = 'regular'
if (ro == 'o'):
    reof = 'offset'

print('Do you want a supercell in the z-direction (2 units)? (enter y/n)')

sup = str(input())

if (sup == 'y'):
    super = 'super'

if (sup == 'n'):
    super = 'unit'

print('Enter additive factor for lengthening the cell in z-direction [0.0]:')

cell_z_factor = input()

if (len(cell_z_factor) == 0):
    cell_z_factor = float('0.0')

cell_z_factor = float(cell_z_factor)

print('Enter inorganic layer separation in units of molecule length [2.0]:')

delta_z_factor = input()

if (len(delta_z_factor) == 0):
    delta_z_factor = float('2.0')

delta_z_factor = float(delta_z_factor)

frame = read('/home/essil/Documents/structure_maker/inorganic/' + cell_type + 'n' + str(n) + '_' + reof + '_' + super + '.traj')

mol = read(str(sys.argv[1]))

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

##### Find positions of nitrogens/carbons in the frame

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

Long = np.split(Long, Long_counter)

NC_vec = np.empty((len(Long), 3))
N_com = np.empty((len(Long), 3))
N_pos = np.empty((len(Long), 3))
counter = 0

for long in Long:
    long_symbols = np.asarray(organic[long].get_chemical_symbols(), dtype='str')
    N_index = np.asscalar(np.flatnonzero(long_symbols == 'N'))
    N_atom = organic[long[N_index]]
    C_indices = np.flatnonzero(long_symbols == 'C')
    N_neighbours = neighborList.get_neighbors(long[N_index])[0]
    
    N_neighbours_symbols = np.asarray(organic[N_neighbours].get_chemical_symbols(), dtype='str')
    N_neighbours_noC = N_neighbours[np.where(N_neighbours_symbols != 'C')]
    
    C_neighbouring = organic[N_neighbours[np.where(N_neighbours_symbols == 'C')][0]]
    
    N_position = N_atom.position
    C_position = C_neighbouring.position
    
    N_pos[counter] = N_position
    NC_vec[counter] = C_position - N_position
    N_com[counter] = organic[long].get_center_of_mass() - N_position
    counter += 1

NC_vec = normalized(NC_vec)
N_com = normalized(N_com)

###### Get molecule N-COM vector

cutOff_mol = neighborlist.natural_cutoffs(mol)
neighborList_mol = neighborlist.NeighborList(cutOff_mol, self_interaction=False, bothways=True)
neighborList_mol.update(mol)

matrix_mol = neighborList.get_connectivity_matrix()

n_components_mol, component_list_mol = sparse.csgraph.connected_components(matrix_mol)

N_mol_index = mol[mol.symbols == 'N'][0].index
N_mol_pos = mol[N_mol_index].position
N_mol_neighbours = mol[neighborList_mol.get_neighbors(N_mol_index)[0]]

C_mol_neighbouring = N_mol_neighbours[N_mol_neighbours.symbols == 'C'][0]
C_mol_pos = C_mol_neighbouring.position

NC_mol_vec = C_mol_pos - N_mol_pos

NC_mol_vec = normalized(NC_mol_vec)[0]

N_mol_com = mol.get_center_of_mass() - N_mol_pos
N_mol_com = normalized(N_mol_com)[0]

##### Separate inorganic frame into top and bottom parts

inorganic_z = inorganic.get_positions()[:,2]

average_z = np.average(inorganic_z)

inorganic_upper = inorganic[np.where(inorganic_z > average_z)[0]]
inorganic_bottom = inorganic[np.where(inorganic_z <= average_z)[0]]

if(len(inorganic_upper) != len(inorganic_bottom)):
    print('Warning: different number of atoms in inorganic top and bottom parts, is this expected?')

inorganic_upper_z = inorganic_upper.get_positions()[:,2]
inorganic_bottom_z = inorganic_bottom.get_positions()[:,2]

inorganic_z_distance = np.amin(inorganic_upper_z) - np.amax(inorganic_bottom_z)

if (sup == 'n'):
    inorganic_z_distance = 0

molecule_length = get_molecule_length(mol)

delta_z = delta_z_factor * molecule_length - inorganic_z_distance

if (sup == 'y'):
    inorganic_upper.positions[:,2] += delta_z

inorganic_all = inorganic_bottom + inorganic_upper

##### Check which MA's and N's belong to the upper and bottom parts

if(n > 1):
    for i in range(len(MA)):
        
        MA_comz = organic[MA[i]].get_center_of_mass()[2]
        
        if(MA_comz < average_z):
            bottom_MA += sort(organic[MA[i]])
        if(MA_comz > average_z):
            upper_MA += sort(organic[MA[i]])
    if (sup == 'y'):
        upper_MA.positions[:,2] += delta_z


if(len(upper_MA) != len(bottom_MA)):
    print('Warning: different number of atoms in bottom and upper MA. Is this expected?')

##### Get rotation matrices

original = mol.get_positions()

new_mol_positions = np.empty((len(N_com), len(mol), 3))

for i in range(len(NC_vec)):
    
    R = get_rotation_matrix(N_mol_com, N_com[i])
    
    for j in range(len(mol)):
        new_mol_positions[i][j] = np.matmul(R, mol.get_positions()[j])

##### Write

if (cell_type == '1'):
    top_two = np.argpartition(N_pos[:,2], -2)[-2:] #move top two N positions first
    for i in top_two:
        if (sup == 'y'):
            N_pos[i][2] += delta_z

if (cell_type == '2'):
    top_four = np.argpartition(N_pos[:,2], -4)[-4:] #move top four N positions first
    for i in top_four:
        if (sup == 'y'):
            N_pos[i][2] += delta_z

for i in range(len(N_pos)):
    
    if(N_pos[i][2] > average_z):
        if (sup == 'y'):
            N_pos[i][2] += delta_z
    
    mol.set_positions(new_mol_positions[i])
    translate = N_pos[i] - mol[N_mol_index].position
    mol.positions += translate
    ordered += mol
    mol.set_positions(original)

if(n > 1):
    
    for i in range(len(bottom_MA)):
        ordered += bottom_MA[i]
    
    for i in range(len(upper_MA)):
        ordered += upper_MA[i]

ordered += inorganic_all[inorganic_all.symbols == 'Pb']
ordered += inorganic_all[inorganic_all.symbols == 'Br']

##### Translate everything so that an atom is in origin and write approximate cell

#min_distance_index = np.argmin(np.linalg.norm(ordered.get_positions(), axis=1))
#translate_all = ordered.get_positions()[min_distance_index]
#ordered.set_positions( ordered.get_positions() - translate_all )

#cell_x = np.amax(ordered.get_positions()[:,0]) + cell_x_factor
#cell_y = np.amax(ordered.get_positions()[:,1]) + cell_y_factor

# Inorganic layers distances must be consistent

bottom_inorganic_z = np.amin(ordered[(ordered.symbols == 'Pb') | (ordered.symbols == 'Br')].get_positions()[:,2])
top_inorganic_z = np.amax(ordered[(ordered.symbols == 'Pb') | (ordered.symbols == 'Br')].get_positions()[:,2])

cell_z = delta_z_factor * molecule_length + top_inorganic_z - bottom_inorganic_z + cell_z_factor

#ordered.set_cell([cell_x, cell_y, cell_z])

ordered.set_cell(frame.get_cell())
ordered.cell[2][2] = cell_z

#cellx = ordered.get_cell()

prefix = str(sys.argv[1]).split(".")[0]

print('Enter prefix of the name of output file [' + cell_type + 'n' + str(n) + prefix + ro +'_' + super + ']')
name = input()

if(len(name) == 0):
    name = cell_type + 'n'  + str(n) + prefix + ro + '_' + super

write(name + '.traj', ordered)
write(name + '.xyz', ordered)

#np.savetxt('cell' + str(n) + '.dat', cellx)
