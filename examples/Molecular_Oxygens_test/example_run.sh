#!/bin/bash

cat << EOF

 An example of a gas phase docked molecule which contains two O atoms that are not part of the zeolite crystal. 
 The script automatically identifies these O atoms and includes them from the mode calculation.
 The script will sort the atoms in the POSCAR file and if a POTCAR file is present it will sort it accordingly.
 The script will output a MODECAR file containing the normal modes of the molecule sorted as the POSCAR and POTCAR.
 
EOF

python ../../scripts/run.py
