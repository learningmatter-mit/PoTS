#!/usr/bin/env python3

import argparse
import sys
import glob
import os
import numpy as np
from pymatgen.core import Structure, Lattice, Element
from utils import DimerUtils


class DimerParser:
    def parse_arguments():
        """
        Parses command-line arguments for the DIMER mode generation script.
        Returns:
            argparse.Namespace: The parsed arguments.
        """
        parser = argparse.ArgumentParser(description="Generate modes file for DIMER calculations.")
        parser.add_argument(
            "--cp2k",
            action="store_true",
            help="Generate the &MODE_VECS part for the CP2K dimer input file.",
        )
        parser.add_argument(
            "--gulp",
            action="store_true",
            help="Generate the vector part for the GULP dimer input file.",
        )
        parser.add_argument(
            "--ase",
            action="store_true",
            help="Generate the mode vector part for the ASE dimer input file.",
        )
        parser.add_argument(
            "--lammps",
            action="store_true",
            help="Generate the displace_atoms part for the LAMMPS dimer input file.",
        )
        return parser.parse_args()

    def parse_xyz_and_vibdisps(xyz_filepath, vibdisps_filepath):
        """
        Parses atomic coordinates from an XYZ file and vibrational displacement vectors from a related file.

        Args:
            xyz_filepath (str): Path to the XYZ file.
            vibdisps_filepath (str): Path to the vibrational displacement file.

        Returns:
            list: List of dictionaries, each containing 'element', 'index', 'coords', and 'vibdisps'.

        Raises:
            ValueError: If any file is improperly formatted or if there is a mismatch between the number of atoms and displacements.
        """

        # Parse XYZ file
        with open(xyz_filepath, "r") as file:
            lines = file.readlines()
            if len(lines) < 2:
                raise ValueError("The XYZ file is too short to contain valid data.")

            gas_data = []
            for idx, line in enumerate(lines[2:], start=0):
                parts = line.split()
                if len(parts) < 4:
                    raise ValueError(f"Invalid line format in XYZ file: {line}")
                element = parts[0]
                try:
                    coord = [float(parts[1]), float(parts[2]), float(parts[3])]
                except ValueError as e:
                    raise ValueError(f"Non-numeric coordinate found in XYZ file: {line}") from e

                gas_data.append(
                    {
                        "element": element,
                        "index": idx,
                        "coords": coord,
                        "O_docked": False,
                    }
                )
                if element == "O":
                    gas_data[-1]["O_docked"] = True

        # Parse vibrational displacements file
        with open(vibdisps_filepath, "r") as file:
            vibdisps = []
            for line in file:
                if line.strip():
                    vibdisps.append([float(x) for x in line.split()])

            if not vibdisps:
                raise ValueError("Vibrational displacements file is empty or improperly formatted.")
            elif len(vibdisps) != len(gas_data):
                raise ValueError("Mismatch between the number of atoms in the XYZ file and the number of vibrational displacements.")

        # Combine XYZ coordinates and vibrational displacements
        for i in range(len(gas_data)):
            gas_data[i]["vibdisps"] = vibdisps[i]

        # Sort the gas data by element
        gas_data = sorted(gas_data, key=lambda x: x["element"])

        return gas_data

    @staticmethod
    def parse_poscar(poscar_path):
        """
        Parses lattice vectors and atomic coordinates from a POSCAR file.

        Args:
            poscar_path (str): Path to the POSCAR file.
            include_molecular_oxygens (list): Indices of oxygen atoms in the docked molecule.
            exclude_catalyst_hydrogens (list): Indices of hydrogen atoms to be excluded.

        Returns:
            tuple: (lattice_vectors, coords_dict, input_to_modefile_correction)
        """

        with open(poscar_path, "r") as f:
            lines = f.readlines()

        # 1. Extract scaling factor (line 2)
        scaling_factor = float(lines[1].strip())

        # 2. Extract lattice vectors (lines 3-5)
        lattice_vectors = np.array([list(map(float, line.strip().split())) for line in lines[2:5]])

        # 3. Extract element symbols (line 6)
        element_symbols = lines[5].strip().split()

        # 4. Extract number of atoms per element (line 7)
        num_atoms = list(map(int, lines[6].strip().split()))

        # 5. Detect if selective dynamics is used (optional lines 8-9)
        coord_start_line = 8
        if lines[7].strip().startswith(("S", "s")):
            print("Selective dynamics detected.")
            coord_start_line += 1  # Skip selective dynamics flag if present

        # 6. Check if coordinates are Direct or Cartesian (line 8/9)
        coord_type = lines[coord_start_line - 1].strip().lower()
        is_direct = "d" in coord_type  # 'Direct' coordinates if 'd' found, otherwise Cartesian

        # 7. Extract atomic positions (starting from line 9/10 based on selective dynamics)
        atom_positions = []
        for i in range(sum(num_atoms)):
            position = list(map(float, lines[coord_start_line + i].strip().split()[:3]))
            atom_positions.append(position)

        # 8. Build lattice and structure
        lattice = Lattice(scaling_factor * lattice_vectors)

        # Create the structure object
        species = []
        for element, count in zip(element_symbols, num_atoms):
            species.extend([Element(element)] * count)

        structure = Structure(lattice, species, atom_positions, coords_are_cartesian=not is_direct)

        sorted_structure = structure.get_sorted_structure(key=lambda site: site.specie.symbol)

        lattice_vectors = structure.lattice.matrix
        molecular_oxygens, catalyst_hydrogens = DimerUtils.indentify_molecular_oxygens_catalyst_hydrogens(sorted_structure)

        excluded_hydrogens, included_oxygens = [], []
        coords_dict = {str(el): [] for el in sorted_structure.species}

        for i, site in enumerate(sorted_structure.sites):
            element = str(site.species_string)
            coord = list(site.coords)
            if element == "O" and molecular_oxygens:
                if i in molecular_oxygens:
                    included_oxygens.append(coord)
                else:
                    coords_dict[element].append(coord)
            elif element == "H" and catalyst_hydrogens:
                if i in catalyst_hydrogens:
                    excluded_hydrogens.append(coord)
                else:
                    coords_dict[element].append(coord)
            else:
                coords_dict[element].append(coord)

        if included_oxygens:
            coords_dict["O_docked"] = included_oxygens
        if excluded_hydrogens:
            coords_dict["H_acid"] = excluded_hydrogens

        return lattice_vectors, coords_dict, molecular_oxygens, catalyst_hydrogens

    def parse_cif(cif_path):
        """
        Parses the atomic coordinates and lattice vectors from a CIF file and tags specific oxygens
        belonging to the reactant molecule as 'O_docked' and hydrogens belonging to the catalyst framework as "H_acid".

        Args:
            cif_path (str): Path to the CIF file.
            include_oxygens (list of int): List of indices for oxygens in the docked molecule.

        Returns:
            lattice_vectors (list): List of lattice vectors.
            coords_dict (dict): Dictionary containing element types as keys and their corresponding coordinates as values.
        """
        structure = Structure.from_file(cif_path)
        sorted_structure = structure.get_sorted_structure(key=lambda site: site.specie.symbol)
        lattice_vectors = structure.lattice.matrix.tolist()
        molecular_oxygens, catalyst_hydrogens = DimerUtils.indentify_molecular_oxygens_catalyst_hydrogens(sorted_structure)

        if not lattice_vectors:
            raise ValueError("The CIF file does not contain valid lattice vectors or atomic coordinates.")

        # Initialize dictionary to store element coordinates
        coords_dict = {str(el): [] for el in sorted_structure.species}
        included_oxygens, excluded_hydrogens = [], []
        for i, site in enumerate(sorted_structure.sites):
            element = site.specie.symbol
            coord = site.coords.tolist()  # Use site.coords.tolist() for Cartesian coordinates
            if element == "O" and molecular_oxygens:
                if i in molecular_oxygens:
                    included_oxygens.append(coord)
                else:
                    coords_dict[element].append(coord)
            elif element == "H" and catalyst_hydrogens:
                if i in catalyst_hydrogens:
                    excluded_hydrogens.append(coord)
                else:
                    coords_dict[element].append(coord)
            else:
                coords_dict[element].append(coord)

        if included_oxygens:
            coords_dict["O_docked"] = included_oxygens
        if excluded_hydrogens:
            coords_dict["H_acid"] = excluded_hydrogens

        return lattice_vectors, coords_dict, molecular_oxygens, catalyst_hydrogens

    def parse_structure_file(work_dir):
        """
        Parses the crystal structure file.

        This function determines which crystal structure file (POSCAR, CIF) always prioritizing the POSCAR.
        It checks if the relevant file exists and is not empty, then parses it to extract the lattice vectors, atomic coordinates,
        and any required modefile corrections.

        Args:
            work_dir (str): The working directory containing the crystal structure files.

        Returns:
            tuple: A tuple containing:
                - lattice (list): A list of lattice vectors.
                - solid_coords (dict): A dictionary of atomic coordinates grouped by element.
                - input_to_modefile_correction (int): Correction factor for modefile indexing.

        Raises:
            FileNotFoundError: If no valid crystal structure file (POSCAR, CIF) is found
                               in the working directory.
            ValueError: If the specified file is empty or improperly formatted.
        """
        try:
            poscar_filepath = f"{work_dir}/POSCAR"
            check_poscar = os.path.isfile(poscar_filepath)
            if check_poscar:
                if DimerUtils.is_file_empty(poscar_filepath):
                    raise ValueError(f"The POSCAR file at {poscar_filepath} is empty.")
                print(f"Parsing POSCAR file from {poscar_filepath}...")
                return DimerParser.parse_poscar(poscar_filepath)
            else:
                cif_filepath = (glob.glob(os.path.join(work_dir, "*.cif")))[0]
                print(f"Parsing CIF file from {cif_filepath}...")
                if DimerUtils.is_file_empty(cif_filepath):
                    raise ValueError(f"The CIF file at {cif_filepath} is empty.")
                return DimerParser.parse_cif(cif_filepath)
        except FileNotFoundError:
            print("FileNotFoundError: You must have either a POSCAR, CIF file as a crystal input structure in your working directory.")
            sys.exit(1)
