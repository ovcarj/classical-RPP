import sys

from ase.io import read, write

mol = str(sys.argv[1])

string = str.split(mol, '.')

atoms = read(mol)

write(string[0] + '.traj', atoms)
