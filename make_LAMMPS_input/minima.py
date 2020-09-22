#!/usr/bin/python

import numpy as np

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

lammps = LAMMPSlib(lmpcmds=cmds, lammps_header=header, amendments=amendments, log_file=None, keep_alive=True, create_atoms=False, create_box=False, boundary=False)

print('Done with reading LAMMPS file.')
#print(header)
#print(cmds)

##################################################################################

atoms.set_calculator(lammps)
print('Calculator set!')

##################################################################################

energy_i = atoms.get_potential_energy()
print('Starting potential energy: ' + str(energy_i))

#opt = FIRE(atoms, trajectory='PREFIX_opt.traj')

#opt.run(fmax=5e-01)

#ecf = ExpCellFilter(atoms)

#qn = FIRE(ecf)

#cell_traj = Trajectory('PREFIX_copt.traj', 'w', atoms)
#qn.attach(cell_traj)
#qn.run(fmax=0.05)

#Temp = 300
#MaxwellBoltzmannDistribution(atoms, temp=Temp*units.kB)

#npt=NPT(atoms, timestep=0.5*units.fs, temperature=Temp*units.kB, externalstress=0.0,
#ttime=25*units.fs,pfactor=0.06*75.**2*units.fs**2, trajectory='NPT.traj')
#npt.run(1000)

energy_f = atoms.get_potential_energy()
print('Final potential energy: ' + str(energy_f))

#print('Starting global optimization...')

hop = MinimaHopping(atoms, timestep=0.5, Ediff0=1.0, T0=50., optimizer=BFGS, 
minima_threshold=5.0e-3, fmax=3e-02, fmax2=0.1, initial_fmax=5e-02, 
mdmin=900, beta1=1.1, beta2=1.1, beta3=1./1.8, externalstress=0., ttime=25., 
pfactor=0.06*75.**2, k1=3., rt1=0.01, k2=10., rt2=0.0, constrain_bool = True)

hop(totalsteps=100)

