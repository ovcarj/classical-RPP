#!/bin/bash

index=1

file=$(realpath written_data.$7)

line_number=$(($(awk '/Atoms/ {print NR}' ${file}) + 2))

LONG_ORGANIC_NO=$1
LONG_ORGANIC_LEN=$2
SHORT_ORGANIC_NO=$3
SHORT_ORGANIC_LEN=$4
INORGANIC_LEN=$5
AVG_Z=$6

for (( i=0; i < $LONG_ORGANIC_NO; i++ ))
do
    awk -v long_organic=$LONG_ORGANIC_LEN -v lnr=$line_number -v ind=$index ' {if(NR>=lnr && NR<lnr+long_organic) {$2=ind; print} else print}' ${file} > tmp && mv tmp ${file}
    index=$(($index + 1))
    line_number=$(($line_number + $2))
    echo $line_number
done

for (( i=0; i < $SHORT_ORGANIC_NO; i++ ))
do
    awk -v MA=$SHORT_ORGANIC_LEN -v lnr=$line_number -v ind=$index ' {if(NR>=lnr && NR<lnr+MA) {$2=ind; print} else print}' ${file} > tmp && mv tmp ${file}
    index=$(($index + 1))
    line_number=$(($line_number + $4))
done

awk -v inorganic=$INORGANIC_LEN -v lnr=$line_number -v ind=$index -v z=$AVG_Z ' {if(NR>=lnr && NR<lnr+inorganic && $7<z) {$2=ind; print} else print}' ${file} > tmp && mv tmp ${file}
index=$(($index + 1))
awk -v inorganic=$INORGANIC_LEN -v lnr=$line_number -v ind=$index -v z=$AVG_Z ' {if(NR>=lnr && NR<lnr+inorganic && $7>z) {$2=ind; print} else print}' ${file} > tmp && mv tmp ${file}

#awk -v type=$TYPES -v dnr=$DIHEDRAL_NR ' {if(NR>=dnr+2 && NR<=dnr+1+type) {t=$3; $3=$4; $4=t; $5="0.0"; print} else print}' ${file} > tmp

