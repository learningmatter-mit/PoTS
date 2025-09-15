#!/usr/bin/env python
import re
import os
import shutil

# Paths
log_file = "../conformer/confgen.log"  # path to your log
xyz_dir = "../conformer"  # where the xyz files are located
xtb_base_dir = "./xtb_runs"  # where folders for xtb runs will be created

os.makedirs(xtb_base_dir, exist_ok=True)

# Read log and extract cluster indices
cluster_indices = []
with open(log_file, "r") as f:
    for line in f:
        m = re.match(r"Cluster (\d+) energy:", line)
        if m:
            cluster_indices.append(int(m.group(1)))

# Process each cluster
for i in cluster_indices:
    # Original xyz filename
    xyz_file = os.path.join(xyz_dir, f"conf_{i + 1}.xyz")
    if not os.path.exists(xyz_file):
        print(f"Warning: {xyz_file} not found!")
        continue

    # Create folder for xtb run
    conf_folder = os.path.join(xtb_base_dir, f"conf_{i + 1}")
    os.makedirs(conf_folder, exist_ok=True)

    # Copy xyz file
    shutil.copy(xyz_file, os.path.join(conf_folder, "gfn_opt.xyz"))

    # Create a minimal xtb.inp
    inp_file = os.path.join(conf_folder, "xtb.inp")
    with open(inp_file, "w") as f:
        f.write("# XTB input file for optimization\n")
        f.write("$opt\n")
        f.write("method = gfn2\n")
        f.write("charge = 0\n")
        f.write("uhf = 0\n")
        f.write("end\n")

print(f"Prepared {len(cluster_indices)} folders for XTB runs.")
