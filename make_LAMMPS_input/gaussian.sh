#!/bin/bash

file=gaussian_in.com

N_ATOM=$(awk 'NR==1 {print $1}' ${file})

tail -${N_ATOM} ${file} > tmp_atoms

cp gaussian_template.com tmp.com

sed -i "8r tmp_atoms" tmp.com

cp tmp.com gaussian_in.com

rm tmp_atoms tmp.com
