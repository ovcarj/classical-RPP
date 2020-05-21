
import numpy as np
from ase.io import read

R = 3.040 #Pb-Br distance

def pretranslate(MA, final):
    rcm = MA.get_center_of_mass()
    T = final - rcm
    MA.positions = MA.positions + T
    return(MA)

inorganic = read('2n2.traj')
MA = read('MA.traj')

final1 = np.array([6.080, 6.080, 12.160])
final2 = np.array([6.080, 0., 12.160])
final3 = np.array([0., 0., 12.160])
final4 = np.array([0., 6.080, 12.160])

atoms = pretranslate(MA, final1) + inorganic

N_position = MA.positions[np.where(MA.symbols == 'N')]
C_position = MA.positions[np.where(MA.symbols == 'C')]
MA.positions[np.where(MA.symbols == 'C')] = N_position
MA.positions[np.where(MA.symbols == 'N')] = C_position

atoms = pretranslate(MA, final2) + atoms

N_position = MA.positions[np.where(MA.symbols == 'N')]
C_position = MA.positions[np.where(MA.symbols == 'C')]
MA.positions[np.where(MA.symbols == 'C')] = N_position
MA.positions[np.where(MA.symbols == 'N')] = C_position

atoms = pretranslate(MA, final3) + atoms

N_position = MA.positions[np.where(MA.symbols == 'N')]
C_position = MA.positions[np.where(MA.symbols == 'C')]
MA.positions[np.where(MA.symbols == 'C')] = N_position
MA.positions[np.where(MA.symbols == 'N')] = C_position

atoms = pretranslate(MA, final4) + atoms

atoms.set_cell(inorganic.get_cell())

atoms.write('2n2_MA.traj')
