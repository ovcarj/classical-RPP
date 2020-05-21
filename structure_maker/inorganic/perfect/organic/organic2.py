import numpy as np
from ase.io import read

def get_rotation_matrix(a,b):
    scalar = np.dot(a,b)
    cross = np.cross(a,b)
    I = np.identity(3)
    V = np.array(([0., -cross[2], cross[1]],[cross[2], 0., -cross[0]],[-cross[1], cross[0], 0.]))
    VV = np.dot(V,V)
    R = I + V + VV * 1. / (1. + scalar)
    return(R)

def pretranslate(mol, final):
    r_N = mol[np.where(mol.symbols == 'N')[0]][0].position
    T = final - r_N
    mol.positions = mol.positions + T
    return(mol)

Pb_Br = 3.040

delta0 = np.array([4.655 - 5.778, 5.385 - 2.583, 11.489 - 13.486])

bza = read('BZA.cif')
frame = read('1n2_up.traj')

#switch x-y-z

bza_x = np.copy(bza.positions[:,0])
bza_y = np.copy(bza.positions[:,1])
bza_z = np.copy(bza.positions[:,2])

bza.positions[:,0] = bza_z
bza.positions[:,1] = bza_x
bza.positions[:,2] = bza_y

#get rotation axis

vec1 = bza[7].position - bza[12].position
vec2 = bza[8].position - bza[16].position

rot_axis = np.cross(vec1, vec2)

##write first set (from top to bottom the order is 4-1-3-2)

N_vec = bza[np.where(bza.symbols == 'N')[0]][0].position
bza_com = bza.get_center_of_mass()

N_com = bza_com - N_vec

final1 = frame[54].position + delta0

atoms = pretranslate(bza, final1) + frame

final2 = final1 + np.sqrt(2.) * Pb_Br #xy
final2[2] = final1[2] #z

atoms = pretranslate(bza, final2) + atoms

final7 = frame[36].position + delta0
final7[2] += 44.

atoms = pretranslate(bza, final7) + atoms

final8 = final7 + np.sqrt(2.) * Pb_Br
final8[2] = final7[2]

atoms = pretranslate(bza, final8) + atoms

##rotate and write second set

bza.rotate(180, v=rot_axis, center='COM')

final3 = frame[43].position - delta0

atoms = pretranslate(bza, final3) + atoms

final4 = final3 + np.sqrt(2.) * Pb_Br #xy
final4[2] = final3[2] #z

atoms = pretranslate(bza, final4) + atoms

final5 = frame[61].position - delta0

atoms = pretranslate(bza, final5) + atoms

final6 = final5 + np.sqrt(2.) * Pb_Br
final6[2] = final5[2]

atoms = pretranslate(bza, final6) + atoms

##

atoms.set_cell(frame.get_cell())
atoms.cell[2][2] = np.max(atoms.positions[:,2])

atoms.write('1n2_up_final.traj')
