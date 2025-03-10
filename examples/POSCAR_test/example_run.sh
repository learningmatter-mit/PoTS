#!/bin/bash

cat << EOF

 An example of a zeolite crystal containing an organic molecule docked. 
 The input is given in POSCAR format.
 The script will sort the atoms in the POSCAR file and if a POTCAR file is present it will sort it accordingly.
 The script will output a MODECAR file containing the normal modes of the molecule sorted as the POSCAR and POTCAR.

EOF

python ../../scripts/run.py
