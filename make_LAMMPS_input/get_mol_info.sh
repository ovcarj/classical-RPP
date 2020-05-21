#!/bin/bash

file=$(realpath $1/data.$1)

ATOMS_NO=$(grep atoms ${file} | awk '{print $1}')
BONDS_NO=$(grep bonds ${file} | awk '{print $1}')
ANGLES_NO=$(grep angles ${file} | awk '{print $1}')
DIHEDRALS_NO=$(grep dihedrals ${file} | awk '{print $1}')

ATOM_TYPES=$(grep "atom types" $file | awk '{print $1}')
BOND_TYPES=$(grep "bond types" $file | awk '{print $1}')
ANGLE_TYPES=$(grep "angle types" $file | awk '{print $1}')
DIHEDRAL_TYPES=$(grep "dihedral types" $file | awk '{print $1}')

ATOMS_NR=$(awk '/Atoms/ {print NR}' ${file})
MASSES_NR=$(awk '/Masses/ {print NR}' ${file})
BOND_C_NR=$(awk '/Bond Coeffs/ {print NR}' ${file})
ANGLE_C_NR=$(awk '/Angle Coeffs/ {print NR}' ${file})
DIHEDRAL_C_NR=$(awk '/Dihedral Coeffs/ {print NR}' ${file})

echo -e $ATOMS_NO > $1_info.dat
echo -e $BONDS_NO >> $1_info.dat
echo -e $ANGLES_NO >> $1_info.dat
echo -e $DIHEDRALS_NO >> $1_info.dat
echo -e $ATOM_TYPES >> $1_info.dat
echo -e $BOND_TYPES >> $1_info.dat
echo -e $ANGLE_TYPES >> $1_info.dat
echo $DIHEDRAL_TYPES >> $1_info.dat

grep Masses $file -A$(($ATOM_TYPES+1)) | tail -$ATOM_TYPES > tmp_$1_Masses.dat
grep "Bond Coeffs" $file -A$(($BOND_TYPES+1)) | tail -$BOND_TYPES > tmp_$1_BondC.dat
grep "Angle Coeffs" $file -A$(($ANGLE_TYPES+1)) | tail -$ANGLE_TYPES > tmp_$1_AngleC.dat
grep "Dihedral Coeffs" $file -A$(($DIHEDRAL_TYPES+1)) | tail -$DIHEDRAL_TYPES > tmp_$1_DihedralC.dat
grep "Pair Coeffs" $file -A$(($ATOM_TYPES+1)) | tail -$ATOM_TYPES > tmp_$1_PairC.dat

grep Atoms $file -A$(($ATOMS_NO+1)) | tail -$ATOMS_NO > tmp_$1_Atoms.dat
grep Bonds $file -A$(($BONDS_NO+1)) | tail -$BONDS_NO > tmp_$1_Bonds.dat
grep Angles $file -A$(($ANGLES_NO+1)) | tail -$ANGLES_NO > tmp_$1_Angles.dat
grep Dihedrals $file -A$(($DIHEDRALS_NO+1)) | tail -$DIHEDRALS_NO > tmp_$1_Dihedrals.dat
