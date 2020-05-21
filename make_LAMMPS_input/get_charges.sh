#!/bin/bash

file=$(realpath written_data.*)

ATOMS_NO=$(grep atoms ${file} | awk '{print $1}')
ATOMS_NR=$(awk '/Atoms/ {print NR}' ${file})

grep Atoms $file -A$(($ATOMS_NO+1)) | tail -$ATOMS_NO | awk '{print $4}' > charges.dat
