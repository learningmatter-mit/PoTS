#!/usr/bin/env python3
import os
import argparse
import math

def distance(coords1, coords2):
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coords1, coords2)))

def prepare_orca_scan(label, orca_opt_base, orca_scan_base, h_index, c_index, final_distance):
    """
    Retrieve ORCA optimized geometries and create ORCA scan input files.
    """
    conf_dir = os.path.join(orca_opt_base, "conf_1")
    xyz_file = os.path.join(conf_dir, "orca_opt.xyz")

    if not os.path.isfile(xyz_file):
        print(f"Error: {xyz_file} not found!")
        return

    with open(xyz_file) as f:
        lines = [line.strip() for line in f.readlines() if line.strip()][2:]  # skip header

    # Prepare ORCA scan folder
    os.makedirs(orca_scan_base, exist_ok=True)

    geometry = []
    for line in lines:
        parts = line.split()
        atom = parts[0]
        coords = [float(x) for x in parts[1:4]]
        geometry.append((atom, coords))

    c_coords = geometry[c_index][1]
    h_coords = geometry[h_index][1]

    initial_distance = distance(c_coords, h_coords)
    n_images = max(2, int(abs(final_distance - initial_distance) / 0.15) + 1)

    # Write ORCA scan input
    scan_file = os.path.join(orca_scan_base, "orca_scan.inp")
    with open(scan_file, "w") as f:
        f.write("! m062x def2-svp Opt TightSCF TightOpt CHELPG\n\n")
        f.write("%MaxCore 1900\n\n")
        f.write("%pal\n nprocs 8\nend\n\n")
        f.write("%scf\n MaxIter 150\nend\n\n")
        f.write("%geom\n Scan\n\n")
        f.write(f"     B {c_index} {h_index} = {initial_distance:.4f}, {final_distance:.4f}, {n_images}\n\n")
        f.write("end\n\n")
        f.write("* xyz 1 1\n")
        f.write("\n".join(lines) + "\n")
        f.write("*\n")
    print(f"Prepared ORCA scan input for {label} in {orca_scan_base}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare ORCA scan input.")
    parser.add_argument("label", choices=["reactant", "product"])
    parser.add_argument("h_index", type=int, help="Index of moving H atom")
    parser.add_argument("c_index", type=int, help="Index of reacting C atom")
    parser.add_argument("final_distance", type=float, help="Final distance for scan (Å)")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    orca_opt_base = os.path.join(base_dir,"orca_opt", f"{args.label}_orca_opt")
    orca_scan_base = os.path.join(base_dir, "orca_scan")

    prepare_orca_scan(args.label, orca_opt_base, orca_scan_base, args.h_index, args.c_index, args.final_distance)
