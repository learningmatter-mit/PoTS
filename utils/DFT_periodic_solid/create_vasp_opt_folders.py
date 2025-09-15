#!/usr/bin/env python
import os
import shutil
import sys
from pymatgen.core import Structure

if len(sys.argv) < 5:
    print("Usage: create_vasp_folders.py <label> <source_docking_folder> <dest_vasp_folder> [incar_file] [kpoints_file]")
    sys.exit(1)

label = sys.argv[1]  # reactant or product
source_base = sys.argv[2]  # docking folder with conf_1, conf_2 ...
dest_base = sys.argv[3]    # where POSCAR and INCAR will go

incar_file = sys.argv[4] if len(sys.argv) > 4 else "INCAR-opt"
kpoints_file = sys.argv[5] if len(sys.argv) > 5 else "KPOINTS"

os.makedirs(dest_base, exist_ok=True)

for folder in sorted(os.listdir(source_base)):
    folder_path = os.path.join(source_base, folder)
    if os.path.isdir(folder_path) and folder.startswith("conf_"):
        cif_file = os.path.join(folder_path, "docked_pose", "0000.cif")

        # Create destination conf_N folder
        dest_folder = os.path.join(dest_base, folder)
        os.makedirs(dest_folder, exist_ok=True)

        if os.path.isfile(cif_file):
            # Convert CIF → POSCAR
            structure = Structure.from_file(cif_file)
            poscar_file = os.path.join(dest_folder, "POSCAR")
            structure.to(fmt="poscar", filename=poscar_file)
            print(f"Converted {cif_file} → {poscar_file}")
        else:
            print(f"Warning: {cif_file} not found in {folder}")

        # Copy INCAR if exists
        if os.path.isfile(incar_file):
            dst_incar = os.path.join(dest_folder, "INCAR")
            shutil.copy2(incar_file, dst_incar)
            print(f"Copied {incar_file} → {dst_incar}")
        else:
            print(f"Warning: {incar_file} not found, skipping.")


        if os.path.isfile(kpoints_file):
            dst_kpoints = os.path.join(dest_folder, "KPOINTS")
            shutil.copy2(incar_file, dst_kpoints)
            print(f"Copied {kpoints_file} → {dst_kpoints}")
        else:
            print(f"Warning: {kpoints_file} not found, skipping.")
