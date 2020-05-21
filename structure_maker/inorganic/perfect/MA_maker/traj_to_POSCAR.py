from ase.io import read, write
import sys

name = str(sys.argv[1])

atoms = read(name)

write('POSCAR' + name.split('.')[0], format='vasp', images=atoms)
