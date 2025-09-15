#!/usr/bin/env python
import os
import shutil
import sys

if len(sys.argv) != 3:
    print("Usage: python retrive_mins_to_dock.py <label> <orca_base_dir>")
    sys.exit(1)

label = sys.argv[1]  # "reactant" or "product"
orca_base_dir = sys.argv[2]  # path to ORCA optimized geometries

# destination docking folder
docking_base_dir = os.getcwd()  # run script inside reactant_docking or product_docking

# files to copy besides the xyz
extra_files = ["zeolite.cif", "docking_job.sh"]

os.makedirs(docking_base_dir, exist_ok=True)

for folder in sorted(os.listdir(orca_base_dir)):
    folder_path = os.path.join(orca_base_dir, folder)
    if os.path.isdir(folder_path) and folder.startswith("conf_"):
        num = folder.split("_")[1]
        conf_folder = os.path.join(docking_base_dir, f"conf_{num}")
        os.makedirs(conf_folder, exist_ok=True)

        # copy ORCA xyz → docking molecule.xyz
        src_xyz = os.path.join(folder_path, "orca_opt.xyz")
        dst_xyz = os.path.join(conf_folder, "molecule.xyz")
        if os.path.isfile(src_xyz):
            shutil.copy2(src_xyz, dst_xyz)
            print(f"Copied {src_xyz} → {dst_xyz}")
        else:
            print(f"Warning: {src_xyz} not found, skipping.")

        # copy extra files
        for f in extra_files:
            src_file = os.path.join(docking_base_dir, f)
            dst_file = os.path.join(conf_folder, f)
            if os.path.isfile(src_file):
                shutil.copy2(src_file, dst_file)
                print(f"Copied {src_file} → {dst_file}")
            else:
                print(f"Warning: {src_file} not found, skipping.")

