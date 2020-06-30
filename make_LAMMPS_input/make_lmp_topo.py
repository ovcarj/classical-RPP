#!/usr/bin/env python

import sys
import numpy as np
import os

from bash import bash

from ase.io import read, write
from ase.io.xyz import write_xyz

def CalculateOnTeran(molecule_prefix, molecule_start, molecule_end):
    
    MOL = atoms[molecule_start:molecule_end]
    n_mol = len(MOL)
    write_xyz(fileobj='gaussian_in.com', images=MOL)
    
    bash('./gaussian.sh')
    
    print('Calculating ' + molecule_prefix + ' on Isabella...')
    
    bash('./teran.sh ' + molecule_prefix + ' ' + str(n_mol))
    
    print('Isabella calculation done.')
    
    return()

def UpdateDatabase():
    
    database = bash('ls ' + database_dir).value().split()
    
    return database

def MakeTypeDictionary(types):
    bash('grep "Masses" -A' + str(types + 1) + ' written_data.' + prefix + ' | tail -' + str(types) + ' > tmp_TypeDict.dat')
    info_Masses = np.loadtxt('tmp_TypeDict.dat')
    dict_string = '{'
    for i in range(0, types):
        mass = info_Masses[i][1]
        if(13.9 < info_Masses[i][1] < 14.1):
            type = '7'
        if(11.9 < info_Masses[i][1] < 12.1):
            type = '6'
        if(0.9 < info_Masses[i][1] < 1.1):
            type = '1'
        if(206. < info_Masses[i][1] < 208.):
            type = '82'
        if(78. < info_Masses[i][1] < 80.):
            type = '35'
        dict_string += str(i+1) + ': ' + type + ','
    dict_string += '}'
    return(dict_string)

def AverageInorganicZ(atoms, inorganic_indices):
    average_z = np.average(atoms.get_positions()[inorganic_indices][:,2])
    return(average_z)

#copy scripts to current dir

script_dir = os.environ["LAMMPS_SCRIPTS_DIR"]

bash('rsync -av ' + script_dir + ' . --exclude=make_lmp_topo.py')

#update database

database_dir = os.environ["MOLECULE_DATABASE_DIR"]

database = UpdateDatabase()

#enter prefixes

prefix = 'BA'        #large organic molecule
small_prefix = 'MA'  #small organic molecule

# Chemical unit info

N = 1          #number of inorganic layers

N_A = 2        #number of large organic spacers (BA, PEA, ...)
N_B = N-1      #number of smaller organic spacers (MA, FA)
N_Pb = N       #number of lead atoms
N_Br = 3*N+1   #number of bromide atoms

# Large organic molecule info (BA)

N_C = 4
N_H = 12
N_N = 1
N_MOL = N_C + N_H + N_N

# Smaller organic molecule info (MA)

N_c = 1
N_h = 6
N_n = 1
N_mol = N_c + N_h + N_n

# Please format the .xyz file in blocks as following:
#    atoms in molecules always in same order!
# Ordering:
# 1. blocks of large organic spacers
# 2. blocks of small organic spacers (relevant for N > 1)
# 3. the inorganic part at the end in any order

# Read .xyz, make sure that first molecule coordinates are unwrapped for gaussian!

xyz = str(bash('find *xyz'))

atoms = read(xyz, format='xyz')

# Determine how many chemical units are in .xyz file

symbols = np.asarray(atoms.get_chemical_symbols(), dtype='str')
positions = atoms.get_positions()

Pb_indices = np.flatnonzero(symbols == 'Pb')
Pb_number = len(Pb_indices)

N_chem = int(Pb_number/N)    #number of chemical units

# Br indices needed for molecular indexing

Br_indices = np.flatnonzero(symbols == 'Br')
inorganic_indices = np.concatenate((Pb_indices, Br_indices))
N_inorganic = len(inorganic_indices)

print('Found', Pb_number, 'Pb. The .xyz should contain', N_chem*2, 'large organic spacers followed by',
N_chem*N_B, 'smaller organic spacers.')

#Check if the large molecule is already in the database

if prefix in database:
    print('Found ' + prefix + ' in database.')
    bash('cp -r ' + database_dir + prefix + ' .')

else:
    print('Could not find ' + prefix + ' in database. Attempting Gaussian/Amber calculation...')
    CalculateOnTeran(prefix, 0, N_MOL)

#Check if the small molecule is already in the database

if(N > 1):
    
    if small_prefix in database:
        print('Found ' + small_prefix + ' in database.')
        bash('cp -r ' + database_dir + small_prefix + ' .')
    
    else:
        print('Could not find ' + small_prefix + ' in database. Attempting Gaussian/Amber calculation...')
        CalculateOnTeran(small_prefix, N_MOL*2*N_chem, N_MOL*2*N_chem+N_mol)

# Get molecular topology and charges info for big molecule

bash("./get_mol_info.sh "  + prefix)

# Get molecular topology and charges info for small molecule

if(N > 1): bash("./get_mol_info.sh "  + small_prefix)

# Load info into numpy arrays for later printing

info = np.loadtxt(prefix + '_info.dat', dtype=int)

Masses = np.loadtxt('tmp_' + prefix + '_Masses.dat')
BondC = np.loadtxt('tmp_' + prefix + '_BondC.dat')
AngleC = np.loadtxt('tmp_' + prefix + '_AngleC.dat')
DihedralC = np.loadtxt('tmp_' + prefix + '_DihedralC.dat')
PairC = np.loadtxt('tmp_' + prefix + '_PairC.dat')
MOL_Atoms = np.loadtxt('tmp_' + prefix + '_Atoms.dat')

Bonds = np.loadtxt('tmp_' + prefix + '_Bonds.dat', dtype=int)
Angles = np.loadtxt('tmp_' + prefix + '_Angles.dat', dtype=int)
Dihedrals = np.loadtxt('tmp_' + prefix + '_Dihedrals.dat', dtype=int)

if(N > 1):
    
    small_info = np.loadtxt(small_prefix + '_info.dat', dtype=int)
    
    small_Masses = np.loadtxt('tmp_' + small_prefix + '_Masses.dat')
    small_BondC = np.loadtxt('tmp_' + small_prefix + '_BondC.dat')
    small_AngleC = np.loadtxt('tmp_' + small_prefix + '_AngleC.dat')
    small_DihedralC = np.loadtxt('tmp_' + small_prefix + '_DihedralC.dat')
    small_PairC = np.loadtxt('tmp_' + small_prefix + '_PairC.dat')
    small_MOL_Atoms = np.loadtxt('tmp_' + small_prefix + '_Atoms.dat')
    
    small_Bonds = np.loadtxt('tmp_' + small_prefix + '_Bonds.dat', dtype=int)
    small_Angles = np.loadtxt('tmp_' + small_prefix + '_Angles.dat', dtype=int)
    small_Dihedrals = np.loadtxt('tmp_' + small_prefix + '_Dihedrals.dat', dtype=int)

else:
    small_info = np.zeros((8), dtype='int')
    small_Masses = []
    small_BondC = []
    small_AngleC = []
    small_DihedralC = []
    small_PairC = []
    small_MOL_Atoms = []
    
    small_Bonds = []
    small_Angles = []
    small_Dihedrals = []

# Print file

file = open('written_data.' + prefix, 'w')

file.write('LAMMPS data file for ' + prefix + '\n')

# info

file.write(str(len(atoms)) + ' atoms\n')
file.write(str(info[1]*N_A*N_chem + small_info[1]*N_chem*N_B) + ' bonds\n')
file.write(str(info[2]*N_A*N_chem + small_info[2]*N_chem*N_B) + ' angles\n')
file.write(str(info[3]*N_A*N_chem + small_info[3]*N_chem*N_B) + ' dihedrals\n\n')

file.write(str(info[4] + small_info[4] + 2) + ' atom types\n')
file.write(str(info[5] + small_info[5]) + ' bond types\n')
file.write(str(info[6] + small_info[6]) + ' angle types\n')
file.write(str(info[7] + small_info[7]) + ' dihedral types\n\n')
file.close()

# Cell, for now it should be provided in file cell*.dat

#cell_path = str(bash('find cell*dat'))
#cell = np.loadtxt(cell_path)
#atoms.set_cell(cell)
write(filename='tmp_cell_data.dat', images=atoms, units='real', format='lammps-data', force_skew=True)
bash('grep xlo tmp_cell_data.dat -A3 >> written_data.' + prefix)

# Masses

file = open('written_data.' + prefix, 'a')
file.write('\n')

file.write('Masses\n\n')

for i in range(len(Masses)):
    file.write(str(i + 1) + ' ' + str(Masses[i][1]) + '\n')

if(N > 1):
    
    for i in range(len(small_Masses)):
        file.write(str(i + 1 + len(Masses)) + ' ' + str(small_Masses[i][1]) + '\n')

#note: Pb is type len(Masses)+len(masses)+1, Br is type len(Masses)+len(masses)+2

file.write(str(len(Masses) + len(small_Masses) + 1) + ' ' + '207.2\n')
file.write(str(len(Masses) + len(small_Masses) + 2) + ' ' + '79.904\n\n')

#Bond Coeffs

file.write('Bond Coeffs\n\n')

for i in range(len(BondC)):
    file.write(str(i + 1) + ' ' + str(BondC[i][1]) + ' ' + str(BondC[i][2]) + '\n')

if(N > 1):
    
    for i in range(len(small_BondC)):
        file.write(str(i + 1 + len(BondC)) + ' ' + str(small_BondC[i][1]) + ' ' + str(small_BondC[i][2]) + '\n')

file.write('\n')

#Angle Coeffs

file.write('Angle Coeffs\n\n')

for i in range(len(AngleC)):
    file.write(str(i + 1) + ' ' + str(AngleC[i][1]) + ' ' + str(AngleC[i][2]) + '\n')

if(N > 1):
    
    for i in range(len(small_AngleC)):
        file.write(str(i + 1 + len(AngleC)) + ' ' + str(small_AngleC[i][1]) + ' ' + str(small_AngleC[i][2]) + '\n')

file.write('\n')

#Dihedral Coeffs

file.write('Dihedral Coeffs\n\n')

for i in range(len(DihedralC)):
    file.write(str(i + 1) + ' ' + str(DihedralC[i][1]) + ' ' + str(int(DihedralC[i][2])) + ' ' + str(int(DihedralC[i][3])) + ' ' + str(DihedralC[i][4]) + '\n')

if(N > 1):
    
    if len(small_DihedralC.shape) == 1:
        small_DihedralC = [small_DihedralC]
    
    for i in range(len(small_DihedralC)):
        file.write(str(i + 1 + len(DihedralC)) + ' ' + str(small_DihedralC[i][1]) + ' ' + str(int(small_DihedralC[i][2])) + ' ' + str(int(small_DihedralC[i][3])) + ' ' + str(small_DihedralC[i][4]) + '\n')

file.write('\n')

#Atoms

file.write('Atoms\n\n')

for i in range(N_MOL*N_chem*N_A):
    file.write(str(i+1) + ' 0 ' + str(int(MOL_Atoms[i%N_MOL][2])) + ' ' + str(MOL_Atoms[i%N_MOL][3]) + ' ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')

if(N > 1):
    
    for i in range(N_MOL*N_chem*N_A, N_MOL*N_chem*N_A + N_mol*N_chem*N_B):
        file.write(str(i+1) + ' 0 ' + str(int(small_MOL_Atoms[(i - N_MOL*N_A*N_chem) % N_mol][2]) + info[4]) + ' ' + str(small_MOL_Atoms[(i - N_MOL*N_A*N_chem) % N_mol][3]) + ' ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')

for i in range(N_MOL*N_chem*N_A + N_mol*N_chem*N_B, len(atoms)):
    if (symbols[i] == 'Pb'):
        chg = 2.0
        typ = str(len(Masses) + len(small_Masses) + 1)
    if (symbols[i] == 'Br'):
        chg = -1.0
        typ = str(len(Masses) + len(small_Masses) + 2)
    if (symbols[i] != 'Pb' and symbols[i] != 'Br'):
        print('Wrong atom type found where Pb or Br expected, exiting.')
        sys.exit()
    
    file.write(str(i+1) + ' 0 ' + str(typ) + ' ' + str(chg) + ' ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')

file.write('\n')

# Bonds

file.write('Bonds\n\n')

for i in range (len(Bonds)*N_A*N_chem):
    file.write(str(i+1) + ' ' + str(Bonds[i%(len(Bonds))][1]) + ' ' + str(Bonds[i%len(Bonds)][2] + N_MOL*int(i/len(Bonds))) + ' ' + str(Bonds[i%(len(Bonds))][3] + N_MOL*int(i/len(Bonds))) + '\n')

if(N > 1):
    
    for i in range (len(Bonds)*N_A*N_chem, len(Bonds)*N_A*N_chem + len(small_Bonds)*N_B*N_chem):
        file.write(str(i+1) + ' ' + str(small_Bonds[(i - len(Bonds)*N_A*N_chem) % len(small_Bonds)][1] + info[5]) + ' ' + str(small_Bonds[(i - len(Bonds)*N_A*N_chem) % (len(small_Bonds))][2] + N_mol*int((i - len(Bonds)*N_A*N_chem) / len(small_Bonds)) + N_MOL*N_A*N_chem) + ' ' + str(small_Bonds[(i - len(Bonds)*N_A*N_chem) % (len(small_Bonds))][3] + N_mol*int((i - len(Bonds)*N_A*N_chem) / len(small_Bonds)) + N_MOL*N_A*N_chem) + '\n')

file.write('\n')

# Angles

file.write('Angles\n\n')

for i in range (len(Angles)*N_A*N_chem):
    file.write(str(i+1) + ' ' + str(Angles[i%(len(Angles))][1]) + ' ' + str(Angles[i%len(Angles)][2] + N_MOL*int(i/len(Angles))) + ' ' + str(Angles[i%(len(Angles))][3] + N_MOL*int(i/len(Angles))) + ' ' + str(Angles[i%(len(Angles))][4] + N_MOL*int(i/len(Angles))) + '\n')

if(N > 1):
    
    for i in range (len(Angles)*N_A*N_chem, len(Angles)*N_A*N_chem + len(small_Angles)*N_B*N_chem):
        file.write(str(i+1) + ' ' + str(small_Angles[(i - len(Angles)*N_A*N_chem) % len(small_Angles)][1] + info[6]) + ' ' + str(small_Angles[(i - len(Angles)*N_A*N_chem) % (len(small_Angles))][2] + N_mol*int((i - len(Angles)*N_A*N_chem) / len(small_Angles)) + N_MOL*N_A*N_chem) + ' ' + str(small_Angles[(i - len(Angles)*N_A*N_chem) % (len(small_Angles))][3] + N_mol*int((i - len(Angles)*N_A*N_chem) / len(small_Angles)) + N_MOL*N_A*N_chem) + ' ' + str(small_Angles[(i - len(Angles)*N_A*N_chem) % (len(small_Angles))][4] + N_mol*int((i - len(Angles)*N_A*N_chem) / len(small_Angles)) + N_MOL*N_A*N_chem) + '\n')

file.write('\n')

# Dihedrals

file.write('Dihedrals\n\n')

for i in range (len(Dihedrals)*N_A*N_chem):
    file.write(str(i+1) + ' ' + str(Dihedrals[i%(len(Dihedrals))][1]) + ' ' + str(Dihedrals[i%len(Dihedrals)][2] + N_MOL*int(i/len(Dihedrals))) + ' ' + str(Dihedrals[i%(len(Dihedrals))][3] + N_MOL*int(i/len(Dihedrals))) + ' ' + str(Dihedrals[i%(len(Dihedrals))][4] + N_MOL*int(i/len(Dihedrals))) + ' ' + str(Dihedrals[i%(len(Dihedrals))][5] + N_MOL*int(i/len(Dihedrals))) + '\n')

if(N > 1):
    
    for i in range (len(Dihedrals)*N_A*N_chem, len(Dihedrals)*N_A*N_chem + len(small_Dihedrals)*N_B*N_chem):
        file.write(str(i+1) + ' ' + str(small_Dihedrals[(i - len(Dihedrals)*N_A*N_chem) % len(small_Dihedrals)][1] + info[7]) + ' ' + str(small_Dihedrals[(i - len(Dihedrals)*N_A*N_chem) % (len(small_Dihedrals))][2] + N_mol*int((i - len(Dihedrals)*N_A*N_chem) / len(small_Dihedrals)) + N_MOL*N_A*N_chem) + ' ' + str(small_Dihedrals[(i - len(Dihedrals)*N_A*N_chem) % (len(small_Dihedrals))][3] + N_mol*int((i - len(Dihedrals)*N_A*N_chem) / len(small_Dihedrals)) + N_MOL*N_A*N_chem) + ' ' + str(small_Dihedrals[(i - len(Dihedrals)*N_A*N_chem) % (len(small_Dihedrals))][4] + N_mol*int((i - len(Dihedrals)*N_A*N_chem) / len(small_Dihedrals)) + N_MOL*N_A*N_chem) + ' ' + str(small_Dihedrals[(i - len(Dihedrals)*N_A*N_chem) % (len(small_Dihedrals))][5] + N_mol*int((i - len(Dihedrals)*N_A*N_chem) / len(small_Dihedrals)) + N_MOL*N_A*N_chem) + '\n')

file.close()

# Pair Coeffs

if(N > 1):
    
    small_PairC[:,0] += info[4]
    small_Masses[:,0] += info[4] #this is to correct the atom types
    
    all_PairC = np.concatenate((PairC, small_PairC))
    all_Masses = np.concatenate((Masses, small_Masses))

else:
    
    all_PairC = PairC
    all_Masses = Masses

file = open('tmp_' + prefix + 'PbBr_PairC.dat', 'w')

for i in range(0, len(all_PairC)):
    for j in range(0, len(all_PairC)):
        if(i==j):
            file.write('pair_coeff      ' + str(int(i+1)) + ' ' + str(int(j+1)) + ' ' + 'lj/cut/coul/long  ' + str(all_PairC[int(i)][1]) + ' ' + str(all_PairC[int(i)][2]) + '\n')
        if(i<j):
            eps=str(np.sqrt(all_PairC[int(i)][1]*all_PairC[int(j)][1]))
            sigma=str(0.5*(all_PairC[int(i)][2]+all_PairC[int(j)][2]))
            file.write('pair_coeff      ' + str(int(i+1)) + ' ' + str(int(j+1)) +  ' ' + 'lj/cut/coul/long  ' + eps + ' ' + sigma + '\n')

h_all_PairC = all_PairC[np.where(all_Masses[:,1] < 2.)] #get hydrogen PairC to write LJ for inorganic-H interaction

Pb_typ = str(info[4] + small_info[4] + 1)
Br_typ = str(info[4] + small_info[4] + 2)

eps_Pb = 0.012484076
sigma_Pb = 3.45998

eps_Br = 0.167296975
sigma_Br = 3.841614286

for i in range(0, len(h_all_PairC)):
    eps = str(np.sqrt(h_all_PairC[i][1] * eps_Pb))
    sigma = str(0.5*(h_all_PairC[i][2] + sigma_Pb))
    file.write('pair_coeff      ' + str(int(h_all_PairC[i][0])) + ' ' + Pb_typ + ' ' + 'lj/cut/coul/long  ' + eps + ' ' + sigma + '\n')

for i in range(0, len(h_all_PairC)):
        eps = str(np.sqrt(h_all_PairC[i][1] * eps_Br))
        sigma = str(0.5*(h_all_PairC[i][2] + sigma_Br))
        file.write('pair_coeff      ' + str(int(h_all_PairC[i][0])) + ' ' + Br_typ + ' ' + 'lj/cut/coul/long  ' + eps + ' ' + sigma + '\n')

file.write('\n')

buckNP = np.array((32690390.937995, 0.150947, 0.000000))
buckNB = np.array((94836.351975893, 0.3352375, 0.000000))

buckCP = np.array((32690390.937995, 0.150947, 0.000000))
buckCB = np.array((94836.351975893, 0.3352375, 0.000000))

buckPP = np.array((74933300.5606326, 0.123246948356808, 0.000000))
buckPB = np.array((110223.38165565, 0.302100469483568, 0.000000))

buckBB = np.array((24274.90558983, 0.45286103286385, 654.4127155))

C_types = all_Masses[np.where(np.abs(all_Masses[:,1] - 12.01) < 0.5 )][:,0]
N_types = all_Masses[np.where(np.abs(all_Masses[:,1] - 14.0067) < 0.5 )][:,0]

for N_typ in N_types:
     
     N_typ = str(int(N_typ))
     
     file.write('pair_coeff      ' + N_typ + ' ' + Pb_typ + ' '
             + 'buck/coul/long    ' + str(buckNP[0]) + ' ' + str(buckNP[1]) + ' ' + str(buckNP[2]) + '\n')
     
     file.write('pair_coeff      ' + N_typ + ' ' + Br_typ + ' '
             + 'buck/coul/long    ' + str(buckNB[0]) + ' ' + str(buckNB[1]) + ' ' + str(buckNB[2]) + '\n')

for C_typ in C_types:
    
    C_typ = str(int(C_typ))
     
    file.write('pair_coeff      ' + C_typ + ' ' + Pb_typ + ' '
            + 'buck/coul/long    ' + str(buckCP[0]) + ' ' + str(buckCP[1]) + ' ' + str(buckCP[2]) + '\n')
     
    file.write('pair_coeff      ' + C_typ + ' ' + Br_typ + ' '
            + 'buck/coul/long    ' + str(buckCB[0]) + ' ' + str(buckCB[1]) + ' ' + str(buckCB[2]) + '\n')

file.write('pair_coeff      ' + Pb_typ + ' ' + Pb_typ + ' '
        + 'buck/coul/long    ' + str(buckPP[0]) + ' ' + str(buckPP[1]) + ' ' + str(buckPP[2]) + '\n')

file.write('pair_coeff      ' + Pb_typ + ' ' + Br_typ + ' '
        + 'buck/coul/long    ' + str(buckPB[0]) + ' ' + str(buckPB[1]) + ' ' + str(buckPB[2]) + '\n')

file.write('pair_coeff      ' + Br_typ + ' ' + Br_typ + ' '
        + 'buck/coul/long    ' + str(buckBB[0]) + ' ' + str(buckBB[1]) + ' ' + str(buckBB[2]))

file.write('\n')

file.close()

bash('sed -i "s/PREFIX/' + prefix + '/g" in_template.file')
bash('cp in_template.file in_' + prefix + 'PbBr.in')
bash('sed -i "s/PREFIX/' + prefix + '/g" minima.py')
bash('sed -i "s/PREFIX/' + prefix + '/g" opt.py')
bash('sed -i "14r tmp_' + prefix + 'PbBr_PairC.dat' + '" in_' + prefix + 'PbBr.in')
bash('./get_charges.sh ' + prefix)

#Make type dictionary

types = info[4] + small_info[4] + 2
Z_of_type_dict = MakeTypeDictionary(types)

bash('sed -i "s/DICT/' + Z_of_type_dict + '/g" minima.py')
bash('sed -i "s/DICT/' + Z_of_type_dict + '/g" opt.py')

#Molecule indexing

if ('super' in xyz):
    average_z = AverageInorganicZ(atoms, inorganic_indices)
elif ('unit' in xyz):
    average_z = -10000.
else:
    print('Warning: strings "super" or "unit" not found in the name of .xyz file. Assuming supercell in z-direction. This may lead to wrong molecular indexing in LAMMPS file.')
    average_z = AverageInorganicZ(atoms, inorganic_indices)

bash('./molecule_indices.sh ' + str(int(N_chem * 2)) + ' ' + str(int(N_MOL)) + ' ' + str(int(N_chem * N_B)) + ' ' + str(int(N_mol)) + ' ' + str(int(N_inorganic)) + ' ' + str(average_z) + ' ' + prefix)

#print('./molecule_indices.sh ' + str(int(N_chem * 2)) + ' ' + str(int(N_MOL)) + ' ' + str(int(N_chem * N_B)) + ' ' + str(int(N_mol)) + ' ' + str(int(N_inorganic)) + ' ' + str(average_z) + ' ' + prefix)

#Cleanup

bash('rm tmp*')
bash('rm *_info.dat')
bash('rm amber2lammps.sh amber2lammps.py gaussian.sh  gaussian_template.com  get_charges.sh  get_mol_info.sh  in_template.file  teran.sh gaussian_in.com neighbor.py molecule_indices.sh README*')
#bash('rm -r ' + prefix)

#if (N > 1):
#    bash('rm -r ' + small_prefix)

bash('mkdir ready')
bash('cp -t ready written_data.' + prefix + ' in_' + prefix + 'PbBr.in' + ' *xyz cell*dat charges.dat minima.py opt.py' )
bash('rm *dat* in* minima.py opt.py')
