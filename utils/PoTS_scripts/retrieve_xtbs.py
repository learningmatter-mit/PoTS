#!/usr/bin/env python
import os
import re

xtb_dir = "../../xtb/product_xtb/xtb_runs"
orca_base_dir = "./orca_opts"
folder_energies = {}

for folder in sorted(os.listdir(xtb_dir)):
    folder_path = os.path.join(xtb_dir, folder)
    out_file = os.path.join(folder_path, "gfn_opt.out")
    if os.path.isdir(folder_path) and os.path.isfile(out_file):
        with open(out_file, "r") as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    # Match the number before 'Eh'
                    match = re.search(r"TOTAL ENERGY\s+([-+]?\d*\.\d+)", line)
                    if match:
                        folder_energies[folder] = float(match.group(1))
                        break  # we only need the first TOTAL ENERGY

if not folder_energies:
    print("No energies found! Check your gfn_opt.out files.")
    exit(1)

# Step 2: Sort by energy and pick two lowest
sorted_folders = sorted(folder_energies.items(), key=lambda x: x[1])
top2 = sorted_folders[:2]

print("Two most stable XTB calculations:")
for folder, energy in top2:
    print(f"{folder}: {energy} Eh")

# Step 3: Extract final geometry and create ORCA input folders
for i, (folder, energy) in enumerate(top2, start=1):
    folder_path = os.path.join(xtb_dir, folder)
    xyz_in = os.path.join(folder_path, "xtbopt.xyz")

    if not os.path.isfile(xyz_in):
        print(f"xtbopt.xyz not found in {folder}")
        continue

    with open(xyz_in, "r") as f:
        lines = f.readlines()

    # skip first two lines (natoms + energy)
    geometry_lines = [line.strip() for line in lines[2:] if line.strip()]

    # Create ORCA folder
    orca_folder = os.path.join(orca_base_dir, f"conf_{i}")
    os.makedirs(orca_folder, exist_ok=True)

    # Save clean xyz for reference
    xyz_file = os.path.join(orca_folder, "xtb_opt.xyz")
    with open(xyz_file, "w") as f:
        f.write("\n".join(geometry_lines) + "\n")

    # Create ORCA input file
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

print(f"Prepared ORCA folders in {orca_base_dir}")
