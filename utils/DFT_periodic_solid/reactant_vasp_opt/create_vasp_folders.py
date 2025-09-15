#!/usr/bin/env python
import os
import shutil
from pymatgen.core import Structure

# Source folder: where conf_1, conf_2, ... with 0000.cif live
source_base = "../../docking/reactant_docking"

# Destination folder: where conf_1, conf_2, ... with POSCAR and INCAR should be created
dest_base = "."

# INCAR to copy
incar_file = "INCAR-opt"  # must be in the same folder where you run this script, or give full path

os.makedirs(dest_base, exist_ok=True)

for folder in sorted(os.listdir(source_base)):
    folder_path = os.path.join(source_base, folder)
    if os.path.isdir(folder_path) and folder.startswith("conf_"):
        cif_file = os.path.join(folder_path, "docked_pose", "0000.cif")

        # Create destination conf_N folder
        dest_folder = os.path.join(dest_base, folder)
        os.makedirs(dest_folder, exist_ok=True)

        if os.path.isfile(cif_file):
            # Convert CIF → POSCAR in destination
            structure = Structure.from_file(cif_file)
            poscar_file = os.path.join(dest_folder, "POSCAR")
            structure.to(fmt="poscar", filename=poscar_file)
            print(f"Converted {cif_file} → {poscar_file}")
        else:
            print(f"Warning: {cif_file} not found in {folder}")

        # Copy INCAR
        if os.path.isfile(incar_file):
            dst_incar = os.path.join(dest_folder, "INCAR")
            shutil.copy2(incar_file, dst_incar)
            print(f"Copied {incar_file} → {dst_incar}")
        else:
            print(f"Warning: {incar_file} not found, skipping.")

