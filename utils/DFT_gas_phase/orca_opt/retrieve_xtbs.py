#!/usr/bin/env python
import os
import re
import argparse

def process_xtb(label, xtb_dir, orca_base_dir):
    folder_energies = {}

    for folder in sorted(os.listdir(xtb_dir)):
        folder_path = os.path.join(xtb_dir, folder)
        if not os.path.isdir(folder_path) or not folder.startswith("conf_"):
            continue

        out_file = os.path.join(folder_path, "gfn_opt.out")
        if not os.path.isfile(out_file):
            print(f"Warning: {out_file} not found, skipping.")
            continue

        with open(out_file, "r") as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    match = re.search(r"TOTAL ENERGY\s+([-+]?\d*\.\d+)", line)
                    if match:
                        folder_energies[folder] = float(match.group(1))
                        break


    if not folder_energies:
        print(f"No energies found in {xtb_dir}! Check your gfn_opt.out files.")
        return

    sorted_folders = sorted(folder_energies.items(), key=lambda x: x[1])
    top2 = sorted_folders[:2]

    print(f"Two most stable XTB calculations for {label}:")
    for folder, energy in top2:
        print(f"{folder}: {energy} Eh")

    # Create ORCA input folders
    for i, (folder, energy) in enumerate(top2, start=1):
        folder_path = os.path.join(xtb_dir, folder)
        xyz_in = os.path.join(folder_path, "xtbopt.xyz")

        if not os.path.isfile(xyz_in):
            print(f"xtbopt.xyz not found in {folder}")
            continue

        with open(xyz_in, "r") as f:
            lines = f.readlines()

        geometry_lines = [line.strip() for line in lines[2:] if line.strip()]

        orca_folder = os.path.join(orca_base_dir, f"conf_{i}")
        os.makedirs(orca_folder, exist_ok=True)

        xyz_file = os.path.join(orca_folder, "xtb_opt.xyz")
        with open(xyz_file, "w") as f:
            f.write("\n".join(geometry_lines) + "\n")

        orca_inp = os.path.join(orca_folder, "orca_opt.inp")
        with open(orca_inp, "w") as f:
            f.write("""! m062x def2-svp Opt TightSCF TightOpt CHELPG

%MaxCore 1900

%pal
 nprocs 8
end

%scf
 MaxIter 150
end

%geom
 MaxIter 150
end

* xyz 1 1
""")
            f.write("\n".join(geometry_lines) + "\n")
            f.write("*\n")

    print(f"Prepared ORCA folders for {label} in {orca_base_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve XTB outputs and prepare ORCA inputs.")
    parser.add_argument("label", choices=["reactant", "product"], help="Which molecule to process")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))  # utils/DFT_gas_phase/orca_opt
    xtb_dir = os.path.abspath(os.path.join(base_dir, f"../../xtb/{args.label}_xtb/"))
    orca_base_dir = os.path.join(base_dir, f"{args.label}_orca_opt")

    os.makedirs(orca_base_dir, exist_ok=True)
    process_xtb(args.label, xtb_dir, orca_base_dir)

