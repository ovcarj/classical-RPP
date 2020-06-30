#!/usr/bin/python

import numpy as np
from numpy.linalg import norm

from ase import Atom, Atoms
from ase.build import bulk
from ase.calculators.lammpslib import LAMMPSlib
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import *
from ase import units
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
from ase.md.nptberendsen import Inhomogeneous_NPTBerendsen as NPTber
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.data import atomic_masses, atomic_numbers
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.optimize.minimahopping import MinimaHopping
from ase.constraints import ExpCellFilter, UnitCellFilter
from ase.io import read, write
from ase.calculators.lammpslib import is_upper_triangular

charges = np.loadtxt('charges.dat')
Z_dict = DICT

lammps_data = 'written_data.PREFIX'

filepath = 'in_PREFIXPbBr.in'

atoms = read(lammps_data, format='lammps-data', units='real', Z_of_type=Z_dict)
atoms.set_pbc(True)

file = open(filepath, 'r')

line = file.readline()
cnt = 1

while line:
    if(len(line.split()) != 0):
        if(line.split()[0] == 'pair_style'):
            end_header = cnt - 1
        if(line.split()[0] == 'pair_modify'):
            end_cmds = cnt
    line = file.readline()
    cnt += 1

file.close()

header = []
cmds = []
amendments= []

file = open(filepath, 'r')

line = file.readline()
cnt = 1

while line:
    if(cnt <= end_header):
        header.append(line)
    if((cnt > end_header) and (cnt <= end_cmds)):
        cmds.append(line)
    if(cnt > end_cmds):
        amendments.append(line)
    line = file.readline()
    cnt += 1

file.close()

for i in range(0, len(charges)):
    amendments.append('set atom ' + str(int(i+1)) + ' charge ' + str(charges[i]))

lammps = LAMMPSlib(lmpcmds=cmds, lammps_header=header, amendments=amendments, log_file='asela.log', keep_alive=True, create_atoms=False, create_box=False, boundary=False)

print('Done with reading LAMMPS file.')
#print(header)
#print(cmds)

##################################################################################

atoms.set_calculator(lammps)
print('Calculator set!')

##################################################################################

energy_i = atoms.get_potential_energy()
print('Starting potential energy: ' + str(energy_i))

opt = BFGS(atoms, trajectory='PREFIX_opt.traj')

opt.run(fmax=5e-01)

ecf = ExpCellFilter(atoms)

qn = BFGS(ecf)

cell_traj = Trajectory('PREFIX_copt.traj', 'w', atoms)
qn.attach(cell_traj)
qn.run(fmax=0.0005)

#tri_mat, coord_transform = convert_cell2(atoms.get_cell())

#if coord_transform is not None:
#    atoms.set_positions([np.matmul(coord_transform, position) for position in atoms.get_positions()])
#    atoms.set_cell(tri_mat.transpose())

#Temp = 300
#MaxwellBoltzmannDistribution(atoms, temp=Temp*units.kB)

#npt=NPT(atoms, timestep=0.5*units.fs, temperature=Temp*units.kB, externalstress=0.0,
#ttime=25*units.fs,pfactor=0.06*75.**2*units.fs**2, trajectory='NPT.traj')
#npt.run(1000)

ecf2 = ExpCellFilter(atoms)

qn2 = BFGS(ecf)

#cell2_traj = Trajectory('PREFIX_copt2.traj', 'w', atoms)
#qn2.attach(cell2_traj)
#qn2.run(fmax=0.001)

#opt2 = BFGS(atoms, trajectory='PREFIX_opt2.traj')

#opt2.run(fmax=0.0005)

energy_f = atoms.get_potential_energy()
print('Final potential energy: ' + str(energy_f))
