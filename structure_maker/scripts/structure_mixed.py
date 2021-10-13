#usr/bin/python

import numpy as np
import sys
import os

from ase.io import read, write
from ase.io.trajectory import Trajectory, TrajectoryWriter
from ase import neighborlist
from ase import Atoms
from ase.build import sort

from scipy import sparse

# Function get_molecule_length(mol):
####
#### Input:  ASE atoms object (molecule)
#### Output: Approximate length of molecule (maximum distance between N and C atoms)

def get_molecule_length(mol):
    N_index = np.where(mol.symbols == 'N')
    C_indices = np.where(mol.symbols == 'C')
    
    length = np.amax(mol.get_distances(N_index, C_indices))
    
    return length

# Function normalized(a, axis=-1, order=2)
####
#### Input:  List of vectors
#### Output: List of normalized (unit) vectors

def normalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

# Function get_rotation_matrix(a,b):
####
#### Input:  Two vectors a and b
#### Output: Rotation matrix that rotates a into b,
#### https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

def get_rotation_matrix(a,b):
    scalar = np.dot(a,b)
    cross = np.cross(a,b)
    I = np.identity(3)
    V = np.array(([0., -cross[2], cross[1]],[cross[2], 0., -cross[0]],[-cross[1], cross[0], 0.]))
    VV = np.dot(V,V)
    
    R = I + V + VV * 1. / (1. + scalar)
    
    return(R)

# First we create empty ASE Atoms objects in which we will write atoms
# 1) ordered ===   This will be our final structure
# 2) upper_MA ===  If we have a supercell, we will write MA molecules from the upper layer here
# 3) bottom_MA === MA from bottom layer. If no supercell, this will be all MA molecules

ordered = Atoms()
upper_MA = Atoms()
bottom_MA = Atoms()

# Factors for lenghtening cell in x and y directions, not used?

cell_x_factor = 0.0
cell_y_factor = 0.0

# We ask the user which unit cell type he wants. The x and y lengths of the cell are measured in the units of Pb-I distance
### cell_type = 1 === cell length is sqrt(2)
### cell_type = 2 === cell length is 2
### cell_type = 3 === cell length is 1

print('Enter unit cell type (1 or 2):')
print('X and Y cell lengths:')
print('1 = sqrt(2) x sqrt(2)')
print('2 = 2x2')
print('3 = 1x1')

cell_type = str(input())

# We ask the users how many inorganic layers he wants.

print('Enter number of inorganic layers (n):')

n = int(input())

# We ask the user if he wants the upper layer to be directly above or offset
# for half xy cell length with regards to the bottom layer.
# This is only relevant if a supercell is picked in the next step.

print('Do you want the cell to be regular or offset? (enter r/o)')

ro = str(input())

if (ro == 'r'):
    reof = 'regular'
if (ro == 'o'):
    reof = 'offset'

# We ask the user if he wants a supercell.

print('Do you want a supercell in the z-direction (2 units)? (enter y/n)')

sup = str(input())

if (sup == 'y'):
    super = 'super'

if (sup == 'n'):
    super = 'unit'

# Here the cell can be lengthened in the z-direction. Otherwise the length will be
# just the difference between the z coordinates of the lowest and highest atoms.

print('Enter additive factor for lengthening the cell in z-direction [0.0]:')

cell_z_factor = input()

if (len(cell_z_factor) == 0):
    cell_z_factor = float('0.0')

cell_z_factor = float(cell_z_factor)

# If a supercell is picked, delta_z_factor*molecule_length will be the distance between top and bottom
# parts of the structure. If the unit cell is chosen, this will be the seperation
# between bottom 2 and upper 2 long organic molecules.

print('Enter inorganic layer separation in units of molecule length [2.0]:')

delta_z_factor = input()

if (len(delta_z_factor) == 0):
    delta_z_factor = float('2.0')

delta_z_factor = float(delta_z_factor)

# This is a code for mixed 50-50 I/Br structures... Should be merged with iodide/bromide.

mix = '_mix'

print('Enter mixed variant (1, 2, 3):')

mix_variant = input()

# Depending on the chosen options,
# template structure is read from the INORGANIC_FRAME_DIR directory

name = cell_type + 'n' + str(n) + '_' + reof + '_' + super + mix + mix_variant

frame = read(os.environ["INORGANIC_FRAME_DIR"] + name + '.traj')

# The corresponding N position / COM pairs are read from .npy files

N_pos = np.load(os.environ["INORGANIC_FRAME_DIR"] + '/N_COM_pairs/' + name + '_N_pos.npy')
com = np.load(os.environ["INORGANIC_FRAME_DIR"] + '/N_COM_pairs/' + name + '_com.npy')

N_com = normalized(com - N_pos)

# The long organic molecule is read from the script argument (e.g. PEA.traj, BZA.traj, ...)

mol = read(str(sys.argv[1]))

# We create lists of indices of different parts of the template structre

symbols = np.asarray(frame.get_chemical_symbols(), dtype='str')

Pb_indices = np.flatnonzero(symbols == 'Pb')
I_indices = np.flatnonzero(symbols == 'I')
Br_indices = np.flatnonzero(symbols == 'Br')
N_indices = np.flatnonzero(symbols == 'N')
C_indices = np.flatnonzero(symbols == 'C')
H_indices = np.flatnonzero(symbols == 'H')

organic_indices = np.concatenate((N_indices, C_indices, H_indices))

inorganic_indices = np.concatenate((Pb_indices, Br_indices, I_indices))

# We create seperate structures containing only organic and inorganic atoms

organic = frame[organic_indices]
inorganic = frame[inorganic_indices]

##### see e.g. https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html#ase.neighborlist.get_connectivity_matrix

if (n>1):
    organic = frame[organic_indices]
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
    
    MA = np.split(MA, MA_counter)

###### Get input molecule (PEA.traj, BZA.traj., ...) N-COM vector

N_mol_index = mol[mol.symbols == 'N'][0].index
N_mol_pos = mol[N_mol_index].position
N_mol_com = mol.get_center_of_mass() - N_mol_pos
N_mol_com = normalized(N_mol_com)[0]

###### Now we separate the inorganic frame from the template into top and bottom parts
###### so we can easily increase/decrease the distance between them to accomodate the
###### long molecules. In the case of unit cell (sup == 'n'), this seperation doesn't
###### do anything because we have only one inorganic layer in the template.

inorganic_z = inorganic.get_positions()[:,2]

average_z = np.average(inorganic_z)

inorganic_upper = inorganic[np.where(inorganic_z > average_z)[0]]
inorganic_bottom = inorganic[np.where(inorganic_z <= average_z)[0]]

# The following warning can be ignored in the case of unit cell.

if(len(inorganic_upper) != len(inorganic_bottom)):
    print('Warning: different number of atoms in inorganic top and bottom parts, is this expected?')

inorganic_upper_z = inorganic_upper.get_positions()[:,2]
inorganic_bottom_z = inorganic_bottom.get_positions()[:,2]

# inorganic_z_distance is the distance between highest inorganic atom in bottom layer
# and lowest inorganic atom in upper layer.

inorganic_z_distance = np.amin(inorganic_upper_z) - np.amax(inorganic_bottom_z)

# If we have unit cell, there will be no upper layer.

if (sup == 'n'):
    inorganic_z_distance = 0

# Approximate molecule length is obtained by the get_molecule_length(mol) function.

molecule_length = get_molecule_length(mol)

# delta_z tells us how much will we have to adjust the inorganic frame distance.

delta_z = delta_z_factor * molecule_length - inorganic_z_distance

# In the case of supercell, the upper layer positions are adjusted by delta_z

if (sup == 'y'):
    inorganic_upper.positions[:,2] += delta_z

# Finally the bottom and the adjusted upper inorganic frames are brought back together.
# In the case of unit cell, we didn't do anything.

inorganic_all = inorganic_bottom + inorganic_upper

##### Since we moved the upper inorganic frame, now we check which MA's from the template
##### belong to the upper and bottom parts so we can adjust their positions accordingly.

if(n > 1):
    for i in range(len(MA)):
        
        MA_comz = organic[MA[i]].get_center_of_mass()[2]
        
        if(MA_comz < average_z):
            bottom_MA += sort(organic[MA[i]])
        if(MA_comz > average_z):
            upper_MA += sort(organic[MA[i]])
    if (sup == 'y'):
        upper_MA.positions[:,2] += delta_z

# The following warning can be ignored in the case of unit cell.

if(len(upper_MA) != len(bottom_MA)):
    print('Warning: different number of atoms in bottom and upper MA. Is this expected?')

##### Now we find rotation matrices between N_com and N_mol_com vectors and apply them to
##### mol positions. We get the list "new_mol_positions" which is a list of rotated
##### mol positions. This way, we get a set of positions for the input molecule
##### (PEA.traj, BZA.traj, ...) such that the reoriented N-COM vectors are the same
##### as the ones in the template structure. Later on, we will translate these
##### "new_mol_positions" so that the positions of the N atoms will be the same as
##### in the template structure.

#we store the original positions for later
original = mol.get_positions()

#empty list for the rotated positions
new_mol_positions = np.empty((len(N_com), len(mol), 3))

for i in range(len(N_com)):
    if (N_com[i][2]) > 0:
        a = np.array([0.0, 0.0, 1.0])
    else:
        a = np.array([0.0, 0.0, -1.0])

    R = get_rotation_matrix(N_mol_com, a)

#for i in range(len(N_com)):
#    
#    R = get_rotation_matrix(N_mol_com, N_com[i])
#    
    for j in range(len(mol)):
        new_mol_positions[i][j] = np.matmul(R, mol.get_positions()[j])

##### Writing of the final structure. For later construction of the potential, it is !!!VERY!!! important that
##### the atoms are written in the following order:
##### 1) Atoms of long molecule 1, followed by atoms of long molecule 2, etc. The atoms of respective long molecules must always be in the same order for every long molecule!
##### 2) Atoms of MA molecule 1, followed by atoms of MA molecule 2, etc.. The atoms of respective MA molecules must always be in the same order for every MA molecule!
##### 3) Pb atoms
##### 4) Br atoms

# We already moved MA's by delta_z, now we have to do it for the N atoms of the long molecules.

if (cell_type == '1'):
    # In the case of supercell, move top two N positions of the long molecules first because they have to move 2*delta_z.
    top_two = np.argpartition(N_pos[:,2], -2)[-2:]
    for i in top_two:
        if (sup == 'y'):
            N_pos[i][2] += delta_z

if (cell_type == '2'):
    # In the case of supercell, move top four N positions of the long molecules first, since they have to move 2*delta_z.
    top_four = np.argpartition(N_pos[:,2], -4)[-4:]
    for i in top_four:
        N_pos[i][2] += delta_z

for i in range(len(N_pos)):
    # Now we move the rest of N atoms (all but the lowest ones are moved here).
    if(N_pos[i][2] > average_z):
        if (sup == 'y'):
            N_pos[i][2] += delta_z
    #### In the following lines, we write the input molecules to their new positions.
    mol.set_positions(new_mol_positions[i])
    # Here the translation is done so the N positions match.
    translate = N_pos[i] - mol[N_mol_index].position
    mol.positions += translate
    # The long molecule is added to the final structure.
    ordered += mol
    # The position of molecule is reset (actually not needed, but not wrong either).
    mol.set_positions(original)

# If n>1, MA's are written now.

if(n > 1):
    
    for i in range(len(bottom_MA)):
        ordered += bottom_MA[i]
    
    for i in range(len(upper_MA)):
        ordered += upper_MA[i]

# The inorganic atoms are written.

ordered += inorganic_all[inorganic_all.symbols == 'Pb']
ordered += inorganic_all[inorganic_all.symbols == 'Br']
ordered += inorganic_all[inorganic_all.symbols == 'I']

##### Finally we have to write the approximate cell.
##### Inorganic layers distances must be consistent, so that when the cell is doubled,
##### the distances between inorganic layers are the same.

bottom_inorganic_z = np.amin(ordered[(ordered.symbols == 'Pb') | (ordered.symbols == 'I') | (ordered.symbols == 'Br')].get_positions()[:,2])
top_inorganic_z = np.amax(ordered[(ordered.symbols == 'Pb') | (ordered.symbols == 'I') | (ordered.symbols == 'Br')].get_positions()[:,2])
top_organic_z = np.amax(ordered[(ordered.symbols == 'N') | (ordered.symbols == 'C') | (ordered.symbols == 'H')].get_positions()[:,2])

if (sup =='y'):
    cell_z = delta_z_factor * molecule_length + top_inorganic_z - bottom_inorganic_z + cell_z_factor

if (sup =='n'):
    cell_z = top_organic_z - bottom_inorganic_z + cell_z_factor

# The x and y cell lengths are taken from the template.

ordered.set_cell(frame.get_cell())

# The z cell length from the template is overwritten with our cell_z.

ordered.cell[2][2] = cell_z

# We save the structure in .xyz and .traj format.

prefix = str(sys.argv[1]).split(".")[0]

print('Enter prefix of the name of output file [' + cell_type + 'n' + str(n) + prefix + ro +'_' + super + mix + mix_variant + ']')
name = input()

if(len(name) == 0):
    name = cell_type + 'n'  + str(n) + prefix + ro + '_' + super + mix + mix_variant

write(name + '.traj', ordered)
write(name + '.xyz', ordered)
