#!/usr/bin/env python3
import os
import sys
import numpy as np

def main():
    if len(sys.argv) != 3:
        print("Usage: python parse_hei_scan.py <orca_scan_folder> <orca_ts_folder>")
        sys.exit(1)

    scan_folder = sys.argv[1]
    ts_folder = sys.argv[2]
    os.makedirs(ts_folder, exist_ok=True)

    dat_file = os.path.join(scan_folder, "orca_opt.relaxscanact.dat")
    if not os.path.isfile(dat_file):
        print(f"ERROR: {dat_file} not found")
        sys.exit(1)

    # Parse energies
    energies = []
    with open(dat_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                _, energy = parts
                energies.append(float(energy))

    if not energies:
        print("ERROR: No energies found in relaxscanact.dat")
        sys.exit(1)

    # Identify highest-energy image
    max_idx = int(np.argmax(energies)) + 1  # +1 since files are .001, .002, ...
    xyz_file = os.path.join(scan_folder, f"orca_opt.{max_idx:03d}.xyz")

    if not os.path.isfile(xyz_file):
        print(f"ERROR: {xyz_file} not found")
        sys.exit(1)

    # Read geometry
    with open(xyz_file, "r") as f:
        lines = f.readlines()[2:]  # skip natoms + comment line
        coords = [line.strip() for line in lines]

    # Write ORCA EVF input
    inp_file = os.path.join(ts_folder, "orca_evf.inp")
    inp_text = """! m062x def2-svp OptTS TightSCF TightOpt CHELPG

%MaxCore 1900

%pal
 nprocs 16
end

%scf
 MaxIter 150
end

%geom
 MaxIter 150
 Calc_Hess true
 Recalc_Hess 5
end

* xyz 1 1
"""
    with open(inp_file, "w") as f:
        f.write(inp_text)
        for c in coords:
            f.write(c + "\n")
        f.write("*\n")

    print(f"ORCA EVF input prepared: {inp_file}")


if __name__ == "__main__":
    main()

