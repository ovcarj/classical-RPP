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
#symbols = np.loadtxt('ordered.xyz', skiprows=2, usecols=0, dtype='str')
Z_dict = DICT

lammps_data = 'written_data.PREFIX'

filepath = 'in_PREFIXPbBr.in'

#xyz = 'BAPbBr_vnl.xyz'
#cell = np.loadtxt('cell_BAPbBr.dat')

atoms = read(lammps_data, format='lammps-data', units='real', Z_of_type=Z_dict)
#atoms.symbols = symbols
masses = atoms.get_masses()
#print(masses)
#print(atoms.get_cell())

#atoms.set_cell(cell)
#atoms.set_initial_charges(charges)
atoms.set_pbc(True)

#print(atoms.get_initial_charges())
#print(atoms.get_cell())

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

#for i in range(0, len(numbers)):
#    amendments.append('mass ' + str(int(i+1)) + ' ' + str(atomic_masses[numbers[i]]))

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

#dyn = FIRE(atoms, trajectory='PREFIX_opt.traj')

#dyn.run(fmax=5e-03)

#ecf = ExpCellFilter(atoms)

#qn = FIRE(ecf)

#traj = Trajectory('PREFIX_copt.traj', 'w', atoms)
#qn.attach(traj)
#qn.run(fmax=0.05)

#atoms2 = read('PREFIX_copt.traj', index=-1)
#atoms2.set_calculator(lammps)
#cell2 = atoms2.get_cell()
#cell2[np.abs(cell2) < 1e-03] = 0.0
#atoms2.set_cell(cell2)

#MaxwellBoltzmannDistribution(atoms, temp=300*units.kB)

#dyn2 = NPTber(atoms, timestep=1 * units.fs, temperature=300, fixcm=True, pressure=1.01325, taut=1e3 * units.fs, taup=1e3 * units.fs, compressibility=0.001)#, taut=1000. * units.fs, pressure=1.01325, taup=1000. * units.fs, compressibility=0.01)

#npt_log = Trajectory('npt_log.traj', 'w', atoms)
#dyn2.attach(npt_log.write)

#dyn2.attach(MDLogger(dyn, atoms, 'md.log', header=False, stress=False, peratom=True, mode="a"), interval=1000)
#dyn2 = NPT(atoms, timestep=1*units.fs, temperature=300, externalstress=1e-5, ttime=1000*units.fs,pfactor=0.05, trajectory='NPT.traj')
#dyn2.run(1000)

#npt_log.close()


#dump_interval = 10000
#traj_file = 'BFGS_log.traj'

#traj_writer = Trajectory(traj_file, 'w', atoms)
#dyn.attach(traj_writer.write, interval=dump_interval)

#dyn.run()

#dyn = VelocityVerlet(atoms, 1 * units.fs, trajectory='tmp.traj')


#dyn.run(100)

energy_f = atoms.get_potential_energy()
print('Final potential energy: ' + str(energy_f))

#print('Starting global optimization...')

hop = MinimaHopping(atoms, timestep=0.5, Ediff0=1.0, T0=50., optimizer=BFGS, minima_threshold=4.0e-3, fmax=1e-02, fmax2=0.1, 
mdmin=200, beta1=1.1, beta2=1.1, beta3=1./1.8, externalstress=0., ttime=25., pfactor=0.06*75.**2, k1=10., rt1=0.01, k2=20., rt2=1.5)
hop(totalsteps=100)

