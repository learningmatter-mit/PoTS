#!/bin/bash

# Runs until the first success is found  Also with Monte Carlo docking
echo ""
echo "This example runs the docking of a Diaryl intermediate into an Al-UTL zeolite framework"
echo ""
echo "If job fails to reach a final pose you can tune --threshold_catan, --threshold and --attempts parameters"
echo ""
echo "Running Monte Carlo docking"
echo ""
python3 ~/voronoi/VOID/scripts/dock.py zeolite.cif molecule.xyz -d mcsuccess -s random -f min_catan_distance -o docked_pose --threshold_catan 3.5 --threshold 1.5 --attempts 200000 
echo "Final pose save to docked_pose folder"
