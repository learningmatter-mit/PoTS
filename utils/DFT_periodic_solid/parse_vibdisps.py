#!/usr/bin/env python3
import os
import sys
import re
import math

def parse_frequencies(lines):
    freq_pattern = re.compile(r"\s*(\d+):\s*([-+]?\d*\.\d+(?:[Ee][-+]?\d+)?)\s+cm")
    freqs = []
    for line in lines:
        m = freq_pattern.match(line)
        if m:
            idx = int(m.group(1))
            val = float(m.group(2))
            freqs.append((idx, val))
    return freqs

_num_re = re.compile(r'[-+]?(?:\d*\.\d+|\d+\.\d*)(?:[Ee][-+]?\d+)?')
def parse_num_token(tok):
    m = _num_re.search(tok)
    if m:
        return float(m.group(0))
    tok2 = tok.strip("(),")
    return float(tok2)

def extract_mode_displacements(lines, target_mode):
    start_idx = next((i for i, l in enumerate(lines) if l.strip().startswith("NORMAL MODES")), None)
    if start_idx is None:
        raise RuntimeError("NORMAL MODES section not found in file.")

    displacements = {}
    current_header = None
    capturing = False
    col_index = None

    i = start_idx + 1
    while i < len(lines):
        line = lines[i].rstrip("\n")
        parts = line.split()
        if not parts:
            i += 1
            continue

        if all(p.isdigit() for p in parts):
            current_header = [int(p) for p in parts]
            if target_mode in current_header:
                capturing = True
                col_index = current_header.index(target_mode)
            else:
                if capturing:
                    break
                capturing = False
                col_index = None
            i += 1
            continue

        if re.match(r'^\s*\d+(:)?', line):
            row_tok = parts[0]
            try:
                row_idx = int(row_tok.rstrip(':'))
            except ValueError:
                i += 1
                continue
            vals = []
            for tok in parts[1:]:
                try:
                    vals.append(parse_num_token(tok))
                except Exception:
                    break
            if capturing and col_index is not None and col_index < len(vals):
                displacements[row_idx] = vals[col_index]
            i += 1
            continue

        if capturing:
            break
        i += 1

    if not displacements:
        raise RuntimeError(f"No displacement data for mode {target_mode} found.")

    max_row = max(displacements.keys())
    n_rows = max_row + 1
    n_atoms = n_rows // 3 if n_rows % 3 == 0 else math.ceil(n_rows / 3)

    atoms = []
    for atom_idx in range(n_atoms):
        x = displacements.get(atom_idx*3, 0.0)
        y = displacements.get(atom_idx*3+1, 0.0)
        z = displacements.get(atom_idx*3+2, 0.0)
        atoms.append((x, y, z))

    return atoms


def main():
    if len(sys.argv) != 3:
        print("Usage: create_vibdisps_from_orca.py <orca_outfile> <dest_vasp_ts_folder>")
        sys.exit(1)

    orca_out = sys.argv[1]
    dest_folder = sys.argv[2]

    if not os.path.isfile(orca_out):
        print(f"ERROR: ORCA output file {orca_out} not found.")
        sys.exit(1)

    if not os.path.isdir(dest_folder):
        print(f"ERROR: Destination folder {dest_folder} not found.")
        sys.exit(1)

    with open(orca_out, "r") as f:
        lines = f.readlines()

    freqs = parse_frequencies(lines)
    if not freqs:
        print("No frequencies found in ORCA output.")
        sys.exit(2)

    lowest_mode, lowest_freq = min(freqs, key=lambda x: x[1])
    if lowest_freq >= 0:
        print(f"WARNING: Lowest frequency mode is not negative ({lowest_freq:.2f} cm^-1).")

    atoms = extract_mode_displacements(lines, lowest_mode)

    vibdisps_path = os.path.join(dest_folder, "vibdisps")
    with open(vibdisps_path, "w") as fo:
        for x, y, z in atoms:
            fo.write(f"{x:12.6f} {y:12.6f} {z:12.6f}\n")

    print(f"vibdisps written to {vibdisps_path}")


if __name__ == "__main__":
    main()

