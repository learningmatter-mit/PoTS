#!/usr/bin/env python3

import os
import glob

from parsers import DimerParser
from modes import DimerModes
from writers import DimerWriter


def main():
    """
    Main function to generate mode files for DIMER calculations.

    This function handles the following steps:
    1. Parses command-line arguments to determine the working directory and input/output files.
    2. Loads gas-phase coordinates from an XYZ file and vibrational displacements from a text file.
    3. Parses the lattice vectors and atomic coordinates from POSCAR or .cif file, always prioritizing POSCAR.
    4. Validates the input files and ensures that the necessary data is available.
    5. Generates solid-state modes based on the parsed input data and vibrational displacements.
    6. Depending on the user-specified output format (VASP, CP2K, GULP, ASE, or LAMMPS), writes the mode data to the appropriate file.

    Raises:
        ValueError: If no valid crystal structure file (POSCAR, CIF) is specified.

    Command-line Arguments:
        --cp2k: Generate the &MODE_VECS part for the CP2K dimer input file.
        --gulp: Generate the vector part for the GULP dimer input file.
        --ase: Generate the mode vector part for the ASE dimer input file.
        --lammps: Generate the displace_atoms part for the LAMMPS dimer input file.
    """
    args = DimerParser.parse_arguments()
    work_dir = os.getcwd()

    # Parse the structure file and load the lattice and atomic coordinates
    lattice, solid_coords, molecular_oxygens, catalyst_hydrogens = DimerParser.parse_structure_file(work_dir)

    # Find any .xyz file in the work_dir
    xyz_files = glob.glob(os.path.join(work_dir, "*.xyz"))
    vibdisps_path = f"{work_dir}/vibdisps"

    if not xyz_files:
        raise FileNotFoundError(f"No .xyz files found in {work_dir}.")
    elif len(xyz_files) > 1:
        raise ValueError(f"Multiple .xyz files found in {work_dir}. Please leave only the one yo want to use in the folder.")
    else:
        gas_coords_path = xyz_files[0]

    if not os.path.exists(vibdisps_path):
        raise FileNotFoundError(f"The file 'vibdisps' was not found in the directory: {work_dir}")

    print(f"Parsing gas coordinates from {gas_coords_path}...")
    print(f"Parsing vibrational displacements from {vibdisps_path}...")
    gas_data = DimerParser.parse_xyz_and_vibdisps(gas_coords_path, vibdisps_path)
    dimer_modes = DimerModes(gas_data, solid_coords, lattice)

    print("Generating solid state modes...")
    modes = dimer_modes.create_modecar(
        args=args,
        molecular_oxygens=molecular_oxygens,
        catalyst_hydrogens=catalyst_hydrogens,
        gas_data=gas_data,
    )

    print("Solid state modes generated.")
    dimer_writer = DimerWriter(modes, solid_coords, lattice)
    dimer_writer.write_output_files(work_dir, args)
