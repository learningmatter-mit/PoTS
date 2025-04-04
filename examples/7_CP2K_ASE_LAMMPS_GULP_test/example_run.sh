#!/bin/bash

cat << EOF

 An example of a zeolite crystal containing an organic molecule docked. 
 The input is given in POSCAR format.
 The script parses the additional arguments in -- and passes them to the run.py script.
 It outputs the input files formatted for CP2K, ASE, LAMMPS, and GULP calculations.

EOF

python ../../scripts/run.py --cp2k --ase --lammps --gulp
