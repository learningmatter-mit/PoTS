#!/bin/bash

cat << EOF

 An example of a zeolite crystal containing acid sites with H atoms that are not part of the docked gas-phase molecule.
 Also, the docked molecule contains two O atoms that are not part of the zeolite crystal.
 The script automatically identifies these acid sites and excludes them from the mode calculation 
 and identifies the O atoms that are not part of the zeolite crystal including them in the mode calculation.
 The script will sort the atoms in the POSCAR file and if a POTCAR file is present it will sort it accordingly.
 The script will output a MODECAR file containing the normal modes of the molecule sorted as the POSCAR and POTCAR.
 
EOF

python ../../scripts/run.py
