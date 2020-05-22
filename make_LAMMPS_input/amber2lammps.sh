#/bin/bash

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh

name=ba-BP86

conda activate python2 && python2 amber2lammps.py ${name} && conda deactivate
