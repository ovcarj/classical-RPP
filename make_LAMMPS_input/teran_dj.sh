#!/bin/bash

ssh jovcar@teran.srce.hr << EOF
	cd /shared/jovcar/gaussian/scripting
	mkdir $1
	cd $1
	mkdir opt
	mkdir charges
EOF

scp gaussian_in.com jovcar@teran.srce.hr:/shared/jovcar/gaussian/scripting/$1/opt

ssh jovcar@teran.srce.hr <<EOF

	cd /shared/jovcar/gaussian/scripting/$1/opt
	cp /shared/jovcar/gaussian/gaussian_submit.sh .
	qsub -sync y gaussian_submit.sh

        cd ../charges
        cp /shared/jovcar/gaussian/ginp-chg_dj.com ginp-chg.com
	cp /shared/jovcar/gaussian/gaussian_submit_chg.sh .
	N2=$(($2+2))
	cat ../opt/gaussian_in.log | tac | grep "Coordinates" -B\${N2} -m1 | tac | tail -${2} | awk '{print \$4, \$5, \$6}' > tmp_coords.dat
	grep "+2 1" ../opt/gaussian_in.com -A${2} | tail -${2} | awk '{print \$1}' > tmp_type.dat
	pr -mts' ' tmp_type.dat tmp_coords.dat > tmp_atoms.dat
	sed -i "8r tmp_atoms.dat" ginp-chg.com
	qsub -sync y gaussian_submit_chg.sh

	module load amber/16
	antechamber -i ginp-chg.log -fi gout -o $1.mol2 -fo mol2 -c resp

	parmchk -i $1.mol2 -f mol2 -o $1.frcmod

	cp /shared/jovcar/gaussian/tleap.in .

	sed -i "s/PRE/$1/g" tleap.in

	tleap -f tleap.in
	rm tmp*

EOF

mkdir $1
cd $1

scp jovcar@teran.srce.hr:/shared/jovcar/gaussian/scripting/$1/charges/\{$1.top,$1.crd,$1.frcmod\} .

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh

conda activate python2 && python2 ../amber2lammps.py ${1} && conda deactivate

file=data.$1

TYPES=$(awk '/dihedral types/ {print $1}' ${file})
DIHEDRAL_NR=$(awk '/Dihedral Coeffs/ {print NR}' ${file})

awk -v type=$TYPES -v dnr=$DIHEDRAL_NR ' {if(NR>=dnr+2 && NR<=dnr+1+type) {t=$3; $3=$4; $4=t; $5="0.0"; print} else print}' ${file} > tmp

mv data.$1 bad_dih_data.$1
mv tmp data.$1
