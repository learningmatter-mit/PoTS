#!/usr/bin/env python
import re
import os
import sys
import shutil

def main():
    if len(sys.argv) != 2:
        print("Usage: retrive_conformers.py <path_to_confgen.log>")
        sys.exit(1)

    log_file = sys.argv[1]
    if not os.path.isfile(log_file):
        print(f"Error: {log_file} not found")
        sys.exit(1)

    # figure out if reactant or product
    if "reactant" in log_file:
        label = "reactant"
    elif "product" in log_file:
        label = "product"
    else:
        print("Error: log path must contain 'reactant' or 'product'")
        sys.exit(1)

    # dirs
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
    conformer_dir = os.path.dirname(log_file)
    xyz_dir = conformer_dir
    xtb_base_dir = os.path.join(base_dir, "xtb",f"{label}_xtb")
    os.makedirs(xtb_base_dir, exist_ok=True)

    # parse log
    cluster_indices = []
    with open(log_file, "r") as f:
        for line in f:
            m = re.match(r"Cluster (\d+) energy:", line)
            if m:
                cluster_indices.append(int(m.group(1)))

    print(f"Found {len(cluster_indices)} clusters for {label}")

    xtb_job_template = os.path.join(base_dir, "xtb", "xtb_job.sh")

    for i in cluster_indices:
        conf_name = f"conf_{i+1}"
        xyz_file = os.path.join(xyz_dir, f"{conf_name}.xyz")
        if not os.path.exists(xyz_file):
            print(f"Warning: {xyz_file} not found!")
            continue

        conf_folder = os.path.join(xtb_base_dir, conf_name)
        os.makedirs(conf_folder, exist_ok=True)

        # copy xyz
        shutil.copy(xyz_file, os.path.join(conf_folder, "gfn_opt.xyz"))

        # copy xtb_job.sh
        if os.path.exists(xtb_job_template):
            shutil.copy(xtb_job_template, conf_folder)
        else:
            print(f"Warning: xtb_job.sh not found at {xtb_job_template}")

        # make xtb.inp
        inp_file = os.path.join(conf_folder, "xtb.inp")
        with open(inp_file, "w") as f:
            f.write("# XTB input file for optimization\n")
            f.write("$opt\n")
            f.write("method = gfn2\n")
            f.write("charge = 0\n")
            f.write("uhf = 0\n")
            f.write("end\n")

    print(f"Prepared {len(cluster_indices)} folders for {label} XTB runs at {xtb_base_dir}")


if __name__ == "__main__":
    main()

