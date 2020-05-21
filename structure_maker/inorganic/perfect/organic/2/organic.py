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

#delta0 = np.array([8.254 - 10.840, 8.067 - 7.938, 5.711 - 8.161])
delta0 = np.array([-8.067 + 7.938, -8.254 + 10.840, 5.711 - 8.161])

pea = read('pea.cif')
frame = read('2n1.traj')

#switch x-y-z

pea_x = np.copy(pea.positions[:,0])
pea_y = np.copy(pea.positions[:,1])
pea_z = np.copy(pea.positions[:,2])

pea.positions[:,0] = -pea_y
pea.positions[:,1] = -pea_x
#pea.positions[:,2] = pea_y

#get rotation axis

vec1 = pea[10].position - pea[15].position
vec2 = pea[17].position - pea[11].position

rot_axis = np.cross(vec1, vec2)

##write first set

N_vec = pea[np.where(pea.symbols == 'N')[0]][0].position
pea_com = pea.get_center_of_mass()

N_com = pea_com - N_vec

final1 = frame[5].position + delta0
final1[2] += 17.

atoms = pretranslate(pea, final1) + frame

pea.euler_rotate(phi=120., center='COM')

final2 = np.copy(final1) #xy
final2[0] -= 2. * Pb_Br #x

atoms = pretranslate(pea, final2) + atoms

pea.euler_rotate(phi=-120., center='COM')

final3 = np.copy(final1) #xy
final3[1] += 2. * Pb_Br #y

atoms = pretranslate(pea, final3) + atoms

pea.euler_rotate(phi=120., center='COM')

final4 = np.copy(final3)
final4[0] -= 2. * Pb_Br

atoms = pretranslate(pea, final4) + atoms

pea.euler_rotate(phi=-120., center='COM')

##rotate and write second set

pea.positions[:,2] = -pea.positions[:,2]

final5 = frame[5].position - delta0

atoms = pretranslate(pea, final5) + atoms

pea.euler_rotate(phi=120., center='COM')

final6 = np.copy(final5) #xy
final6[0] -= 2. * Pb_Br #x

atoms = pretranslate(pea, final6) + atoms

pea.euler_rotate(phi=-120., center='COM')

final7 = np.copy(final5) #xy
final7[1] += 2. * Pb_Br #y

atoms = pretranslate(pea, final7) + atoms

pea.euler_rotate(phi=120., center='COM')

final8 = np.copy(final7)
final8[0] -= 2. * Pb_Br

atoms = pretranslate(pea, final8) + atoms

pea.euler_rotate(phi=-120., center='COM')


##
atoms.positions[:,2] -= 2. * Pb_Br
atoms.set_cell(frame.get_cell())
atoms.cell[0][0] += Pb_Br
atoms.cell[1][1] += Pb_Br
atoms.cell[2][2] = np.max(atoms.positions[:,2])

atoms.write('2n1_final.traj')
