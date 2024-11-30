import numpy as np
from pymatgen.io.vasp import Potcar
import os


class DimerWriter:
    def __init__(self, modes, solid_coords, lattice):
        self.modes = modes
        self.solid_coords = solid_coords
        self.lattice = lattice

    def write_vasp(self, output_dir):
        """
        Writes the MODECAR file for VASP dimer calculations.

        The function generates a MODECAR file that contains the vibrational
        displacement vectors calculated for the system. The file is saved in
        the specified output directory with the filename "MODECAR".

        Args:
            output_dir (str): The directory where the MODECAR file will be saved.

        Raises:
            IOError: If there is an error saving the MODECAR file to the specified output directory.
        """
        modecar_path = f"{output_dir}/MODECAR"

        try:
            np.savetxt(modecar_path, np.array(self.modes), fmt="%.8f")
            print(f"MODECAR file successfully saved. {len(self.modes)} displacement vectors were written to {modecar_path}.")
            atom_types_str, total_atoms = self.write_poscar(self.solid_coords, self.lattice)
            atom_types = " ".join(atom_types_str)
            print(f"POSCAR file successfully saved. {total_atoms} atoms were written sorted by atom type: {atom_types}.")
            if os.path.exists("POTCAR"):
                DimerWriter.read_and_sort_potcar("POTCAR", atom_types_str, "POTCAR")
                print(f"POTCAR file sorted alphabetically with atom types as POSCAR: {atom_types}.")
            else:
                print("No POTCAR file found in the current directory.")
                print("Reminder, your POSCAR has been sorted by atom type. Be aware of the POTCAR order when launching VASP.")

        except IOError as e:
            self.handle_io_error(e, modecar_path, "MODECAR file")

    def write_poscar(self, solid_data, lattice_vectors, coordinate_type="Cartesian"):
        """
        Writes a POSCAR file for VASP input from a sorted list of atom dictionaries.

        Parameters:
        filename (str): The name of the POSCAR file to write.
        comment (str): A comment line for the POSCAR file.
        lattice_vectors (list of list of float): The 3x3 matrix of lattice vectors.
        solid_data (list of dict): List of dictionaries with sorted atom symbol and coordinates.
                                  Example: [{"symbol": "Si", "coords": [x, y, z]}, ...]
        coordinate_type (str): Either "Cartesian" or "Direct", specifying the coordinate system.
        """

        # Sort and handle the 'O_docked' condition
        cleaned_atom_coords_dict = {}

        for symbol, coords in solid_data.items():
            # Replace 'O_docked' with 'O' and 'H_acid' with 'H'
            if symbol == "O_docked":
                symbol = "O"
            elif symbol == "H_acid":
                symbol = "H"

            # Append coordinates to the new dictionary
            if symbol in cleaned_atom_coords_dict:
                cleaned_atom_coords_dict[symbol].extend(coords)
            else:
                cleaned_atom_coords_dict[symbol] = coords

        # Sort the atom symbols
        sorted_symbols = sorted(cleaned_atom_coords_dict.keys())

        # Get atom counts and coordinates
        atom_counts = [len(cleaned_atom_coords_dict[symbol]) for symbol in sorted_symbols]
        atom_coordinates = [coords for symbol in sorted_symbols for coords in cleaned_atom_coords_dict[symbol]]

        # Write the POSCAR file
        with open("POSCAR", "w") as f:
            # Write the comment line
            f.write("POSCAR file created from dimer generator script with atoms sorted \n")
            # Write the scale factor (assumed to be 1.0)
            f.write("1.0\n")
            # Write the lattice vectors
            for vector in lattice_vectors:
                f.write(f"{' '.join(f'{x:.16f}' for x in vector)}\n")
            # Write the atomic symbols and counts
            f.write(f"{' '.join(sorted_symbols)}\n")
            f.write(f"{' '.join(str(count) for count in atom_counts)}\n")
            # Write coordinate type (Direct or Cartesian)
            f.write(f"{coordinate_type}\n")
            # Write atomic coordinates
            for coords in atom_coordinates:
                f.write(f"{' '.join(f'{x:.16f}' for x in coords)}\n")
        return sorted_symbols, sum(atom_counts)

    def read_and_sort_potcar(potcar_filename, sorted_atom_symbols, output_potcar_filename):
        """
        Reads a POTCAR file, splits it into sections based on atomic symbols,
        sorts the sections alphabetically according to the provided sorted atom list,
        and writes the sorted POTCAR to a new file.

        Parameters:
        potcar_filename (str): The path to the POTCAR file.
        sorted_atom_symbols (list of str): The sorted list of atom symbols (e.g., ['Al', 'C', 'H', 'O']).
        output_potcar_filename (str): The path to save the sorted POTCAR file.
        """

        # Read the entire POTCAR file
        with open(potcar_filename, "r") as f:
            potcar_content = f.read()

        # Split the POTCAR file into sections by finding the "End of Dataset" marker
        potcar_sections = potcar_content.split("End of Dataset\n")

        # Remove any empty section that may have been created during the split
        potcar_sections = [section for section in potcar_sections if section.strip()]

        # Create a dictionary to map each atomic symbol to its corresponding POTCAR section
        potcar_dict = {}

        for section in potcar_sections:
            # Find the atomic symbol in the header of each section (e.g., "PAW_PBE Si 08Jan2001")
            for line in section.splitlines():
                if "TITEL" in line.strip():
                    # Extract the atomic symbol from the TITEL line
                    atomic_symbol = line.split()[3]
                    potcar_dict[atomic_symbol] = section + "End of Dataset\n"
                    break

        # Create a new list of POTCAR sections sorted by the provided atomic symbols
        sorted_potcar = [potcar_dict[symbol] for symbol in sorted_atom_symbols if symbol in potcar_dict]
        # Write the sorted sections to a new POTCAR file
        with open(output_potcar_filename, "w") as f:
            f.writelines(sorted_potcar)

    def write_cp2k(self, output_dir):
        """
        Writes the MODE_VECS file for CP2K dimer calculations.

        The function generates a file containing vibrational displacement vectors formatted for
        CP2K input. Each vector is labeled as a "VECTOR" line, with an example displacement comment.
        The file is saved in the specified output directory with the filename "CP2KMODEVECS".

        Args:
            output_dir (str): The directory where the CP2KMODEVECS file will be saved.

        Raises:
            IOError: If there is an error saving the CP2KMODEVECS file to the specified output directory.
        """
        modevecs_path = f"{output_dir}/CP2KMODEVECS"
        print(f"Saving &MODE_VECS info to {modevecs_path}...")
        try:
            with open(modevecs_path, "w") as f:
                for mode in self.modes:
                    mode_line = "VECTOR " + " ".join(f"{val:.8f}" for val in mode)
                    f.write(mode_line + "\n")
            print(f"CP2KMODEVECS file successfully saved. {len(self.modes)} displacement vectors were written to {modevecs_path}.")
        except IOError as e:
            self.handle_io_error(e, modevecs_path, "CP2KMODEVECS file")

    def write_gulp(self, output_dir):
        """
        Writes the vibrational displacement vectors to a GULP input file.

        The function generates a file containing vibrational displacement vectors formatted for
        GULP input. Each line contains an element symbol followed by its corresponding displacement
        vector components. The file is saved in the specified output directory with the filename "GULPvector".

        Args:
            output_dir (str): The directory where the GULPvector file will be saved.

        Raises:
            IOError: If there is an error saving the GULPvector file to the specified output directory.
        """
        gulp_path = f"{output_dir}/GULPvector"
        print(f"Saving vector info to {gulp_path}...")
        try:
            with open(gulp_path, "w") as f:
                index = 0
                for element, coords in self.solid_coords.items():
                    for _ in coords:
                        # print(element, _, index)
                        if element == "O_docked" or element == "H_acid":
                            continue
                        else:
                            mode_line = f"{element:<2} " + " ".join(f"{val:.8f}" for val in self.modes[index])
                            f.write(mode_line + "\n")
                            index += 1
            print(f"GULPvector input file successfully saved. {len(self.modes)} displacement vectors were written to {gulp_path}.")
        except IOError as e:
            self.handle_io_error(e, gulp_path, "GULPvector input file")

    def write_ase(self, output_dir):
        """
        Writes the vibrational displacement vectors to an ASE input file.

        The function generates a file containing vibrational displacement vectors formatted for
        use in ASE (Atomic Simulation Environment). The displacements are written in Python list format,
        with each element followed by its corresponding displacement vector. The file is saved in the
        specified output directory with the filename "ASE_INPUT".

        Args:
            output_dir (str): The directory where the ASE input file will be saved.

        Raises:
            IOError: If there is an error saving the ASE input file to the specified output directory.
        """
        ase_path = f"{output_dir}/ASE_INPUT"
        print(f"Saving ASE input info to {ase_path}...")
        try:
            with open(ase_path, "w") as f:
                f.write("mode = [\n")
                index = 0
                for element, coords in self.solid_coords.items():
                    for i, _ in enumerate(coords):
                        mode_line = f"        ({self.modes[index][0]:.8f}, {self.modes[index][1]:.8f}, {self.modes[index][2]:.8f}),"
                        f.write(mode_line + "\n")
                        index += 1
                f.write("]\n")
            print(f"ASE input file successfully saved. {len(self.modes)} displacement vectors were written to {ase_path}.")
        except IOError as e:
            self.handle_io_error(e, ase_path, "ASE input file")

    def write_lammps(self, output_dir):
        """
        Writes the vibrational displacement vectors to a LAMMPS input file.

        The function generates a file containing vibrational displacement vectors formatted for use
        in LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator). The displacements
        are written as commands to displace atoms in specific directions. Each atom's displacement
        is listed with its index in the structure. The file is saved in the specified output directory
        with the filename "LAMMPS_INPUT".

        Args:
            output_dir (str): The directory where the LAMMPS input file will be saved.

        Raises:
            IOError: If there is an error saving the LAMMPS input file to the specified output directory.
        """
        lammps_path = f"{output_dir}/LAMMPS_INPUT"
        print(f"Saving LAMMPS input info to {lammps_path}...")
        try:
            with open(lammps_path, "w") as f:
                atom_index = 1
                for element, coords in self.solid_coords.items():
                    for i, _ in enumerate(coords):
                        displacement_line = f"displace_atoms {atom_index} move {self.modes[atom_index-1][0]:.8f} {self.modes[atom_index-1][1]:.8f} {self.modes[atom_index-1][2]:.8f}\n"
                        f.write(displacement_line)
                        atom_index += 1
            print(f"LAMMPS input file successfully saved. {len(self.modes)} displacement vectors were written to {lammps_path}.")
        except IOError as e:
            self.handle_io_error(e, lammps_path, "LAMMPS input file")

    def write_output_files(dimer_writer, output_dir, args):
        """
        Writes the appropriate output files for DIMER calculations based on the provided arguments.

        This function checks the command-line arguments to determine which output files should be
        generated by the `DimerWriter` instance. It supports various formats including VASP, CP2K,
        GULP, ASE, and LAMMPS. Depending on the arguments, the corresponding output file is written
        to the specified directory.

        Args:
            dimer_writer (DimerWriter): An instance of the `DimerWriter` class, initialized with
                                        the modes and solid coordinates.
            output_dir (str): The directory where the output files will be saved.
            args (argparse.Namespace): Parsed command-line arguments that indicate which output files
                                       to generate.

        Output Files:
            - MODECAR (for VASP) if `--vasp` is specified.
            - CP2KMODEVECS (for CP2K) if `--cp2k` is specified.
            - GULPvector (for GULP) if `--gulp` is specified.
            - ASE_INPUT (for ASE) if `--ase` is specified.
            - LAMMPS_INPUT (for LAMMPS) if `--lammps` is specified.

        Raises:
            IOError: If there is an issue writing any of the files to the specified output directory.
        """
        if not args.cp2k and not args.gulp and not args.ase and not args.lammps:
            dimer_writer.write_vasp(output_dir)
        if args.cp2k:
            dimer_writer.write_cp2k(output_dir)
        if args.gulp:
            dimer_writer.write_gulp(output_dir)
        if args.ase:
            dimer_writer.write_ase(output_dir)
        if args.lammps:
            dimer_writer.write_lammps(output_dir)

    def handle_io_error(self, exception, file_path, file_description="file"):
        """
        Handles IOError exceptions by raising a more descriptive error message.

        Args:
            exception (IOError): The caught IOError exception.
            file_path (str): The path to the file that caused the error.
            file_description (str): A brief description of the file type (default is "file").

        Raises:
            IOError: A re-raised IOError with an enhanced error message.
        """
        raise IOError(
            f"Failed to save {file_description} to {file_path}. Ensure the directory is writable, the path exists, and there are no conflicts with existing files."
        ) from exception
