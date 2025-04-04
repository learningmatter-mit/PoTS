#!/bin/bash

cat << EOF

 An example of a LaMnO catalyst containing an organic molecule docked. 
 The input is given in POSCAR format.
 The script first changes the CATALYST_FRAMEWORK_ELEMENTS line in scripts/chemistry composition to the necessary atomtypes.
 The script will sort the atoms in the POSCAR file. 
 The POTCAR file is now used in this case because this is just a demo example for non-zeolite catalysis uses.
 The script will output a MODECAR file containing the normal modes of the molecule sorted as the POSCAR.
 Finally, the script returns the CATALYST_FRAMEWORK_ELEMENTS line to its original state for zeolites catalysis.

EOF

CHEMISTRY_COMPOSITION="/home/paufv/clean_notebooks/public_dimer_script/scripts/chemistry_composition.py"
LINE_NUMBER=12
CATALYST_COMPOSITION="CATALYST_FRAMEWORK_ELEMENTS = ['La', 'Mn', 'O']"
ZEOLITE_COMPOSITION="CATALYST_FRAMEWORK_ELEMENTS = ['Al', 'Si', 'O']"

echo "Changing CATALYST_FRAMEWORK_ELEMENTS to ['La', 'Mn', 'O']"
echo ""
echo ""
sed -i "${LINE_NUMBER}s/.*/${CATALYST_COMPOSITION}/" "$CHEMISTRY_COMPOSITION"

python ../../scripts/run.py

echo ""
echo ""
echo "Changing CATALYST_FRAMEWORK_ELEMENTS to ['Al', 'Si', 'O'] again"
sed -i "${LINE_NUMBER}s/.*/${ZEOLITE_COMPOSITION}/" "$CHEMISTRY_COMPOSITION"
