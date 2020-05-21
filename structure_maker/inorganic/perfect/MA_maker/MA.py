
import numpy as np
from ase.io import read

R = 3.040 #Pb-Br distance

def pretranslate(MA, final):
    rcm = MA.get_center_of_mass()
    T = final - rcm
    MA.positions = MA.positions + T
    return(MA)

inorganic = read('1n2.traj')
MA = read('MA.traj')

final1 = R * np.array([np.sqrt(2), 0., 2.])
final2 = R * np.array([0., np.sqrt(2), 2.])

atoms = pretranslate(MA, final1) + inorganic

N_position = MA.positions[np.where(MA.symbols == 'N')]
C_position = MA.positions[np.where(MA.symbols == 'C')]
MA.positions[np.where(MA.symbols == 'C')] = N_position
MA.positions[np.where(MA.symbols == 'N')] = C_position

atoms = pretranslate(MA, final2) + atoms

atoms.set_cell(inorganic.get_cell())

atoms.write('1n1_MA.traj')
