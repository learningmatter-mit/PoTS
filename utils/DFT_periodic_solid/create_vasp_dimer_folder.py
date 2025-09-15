#!/usr/bin/env python
import os
import sys
from pymatgen.core import Structure
import shutil

if len(sys.argv) < 3:
    print("Usage: create_vasp_dimer_folders.py  <source_docking_folder> <dest_vasp_folder> [incar_file] [kpoints_file]")
    sys.exit(1)

source_base = sys.argv[1]  # docking folder with conf_1, conf_2 ...
dest_base = sys.argv[2]    # where POSCARs will go

# Optional INCAR / KPOINTS
incar_file = sys.argv[3] if len(sys.argv) > 3 else None
kpoints_file = sys.argv[4] if len(sys.argv) > 4 else None

os.makedirs(dest_base, exist_ok=True)

folder_path = os.path.join(source_base, "conf_1")
cif_file = os.path.join(folder_path, "docked_pose", "0000.cif")
molecule_xyz_file = os.path.join(folder_path, "molecule.xyz")

if not os.path.isfile(cif_file):
    print(f"Warning: {cif_file} not found in {folder_path}")
    

# Create destination conf_N folder
dest_folder = os.path.join(dest_base)
os.makedirs(dest_folder, exist_ok=True)

# Copy INCAR if provided
if incar_file and os.path.isfile(incar_file):
    dst_incar = os.path.join(dest_folder, "INCAR")
    shutil.copy2(incar_file, dst_incar)
    print(f"Copied {incar_file} → {dst_incar}")

# Copy KPOINTS if provided
if kpoints_file and os.path.isfile(kpoints_file):
    dst_kpoints = os.path.join(dest_folder, "KPOINTS")
    shutil.copy2(kpoints_file, dst_kpoints)
    print(f"Copied {kpoints_file} → {dst_kpoints}")

# Copy .cif file
if cif_file and os.path.isfile(cif_file):
    dst_cif = os.path.join(dest_folder, "ts_docked.cif")
    shutil.copy2(cif_file, dst_cif)
    print(f"Copied {cif_file} → {dst_cif}")

# Copy .xyz file
if molecule_xyz_file and os.path.isfile(molecule_xyz_file):
    dst_xyz = os.path.join(dest_folder, "molecule.xyz")
    shutil.copy2(molecule_xyz_file, dst_xyz)
    print(f"Copied {molecule_xyz_file} → {dst_xyz}")
