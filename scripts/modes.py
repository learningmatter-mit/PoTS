import numpy as np
from chemistry_composition import (
    CATALYST_FRAMEWORK_ELEMENTS,
    CATALYST_NON_FRAMEWORK_ELEMENTS,
)


class DimerModes:
    def __init__(
        self,
        gas_data,
        solid_coords,
        lattice_vectors,
    ):
        self.gas_data = gas_data  # List of gas phase coordinates and vibrationial displacements
        self.solid_coords = solid_coords  # List of solid phase coordinates
        self.lattice_vectors = lattice_vectors  # List of lattice vectors

    def adjust_to_unit_cell(self):
        """
        Adjusts both coordinates and vibrational displacements to be within the unit cell.

        Args:
            coords_dict (dict): Dictionary with atomic coordinates and vibrational displacements.

        Returns:
            dict: Adjusted coordinates and vibrational displacements.

        Raises:
            ValueError: If the gas coordinates dictionary is empty or the lattice vectors are not a 3x3 matrix.
        """
        if not self.gas_data:
            raise ValueError("Gas coordinates dictionary is empty.")

        if len(self.lattice_vectors) != 3 or any(len(v) != 3 for v in self.lattice_vectors):
            raise ValueError("Lattice vectors must be a 3x3 matrix.")

        a, b, c = [np.linalg.norm(v) for v in self.lattice_vectors]

        adjusted_data = []
        for entry in self.gas_data:
            adjusted_data.append(
                {
                    "element": entry["element"],
                    "coords": [
                        entry["coords"][0] % a,
                        entry["coords"][1] % b,
                        entry["coords"][2] % c,
                    ],
                    "vibdisps": [
                        entry["vibdisps"][0] % a,
                        entry["vibdisps"][1] % b,
                        entry["vibdisps"][2] % c,
                    ],
                }
            )
        return adjusted_data

    def get_docked_molecule_coords(self, args):
        """
        Extracts coordinates of non-framework atoms.

        Args:
            parsed_poscar_data (dict): Parsed POSCAR data.

        Returns:
            dict: Non-framework coordinates.
        """
        catalyst_elements = CATALYST_FRAMEWORK_ELEMENTS
        if "H_acid" in self.solid_coords.keys():
            catalyst_elements += ["H_acid"]

        return {element: coords for element, coords in self.solid_coords.items() if element not in catalyst_elements}

    def get_docking_rotations(self, adjusted_gas_coords, dock_gas_coords):
        """
        Generates a rotation matrix to align the gas-phase transition state (TS) coordinates
        with the docked structure's coordinates.

        Args:
            adjusted_gas_coords (dict): A dictionary containing element symbols as keys and
                                        lists of adjusted gas coordinates as values.
            dock_gas_coords (dict): A dictionary containing element symbols as keys and
                                    lists of coordinates for the docked structure as values.

        Returns:
            np.ndarray: A 3x3 rotation matrix that aligns the gas-phase TS with the docked structure.

        Raises:
            ValueError: If one or both coordinate sets are empty.
        """

        # Sort both sets of coordinates by element to ensure consistent ordering
        adjusted_gas_coords = sorted(adjusted_gas_coords, key=lambda x: x["element"])
        dock_gas_coords = {k: v for k, v in sorted(dock_gas_coords.items())}

        # Convert dictionaries to lists of coordinates
        adjusted_gas_coords = np.array([entry["coords"] for entry in adjusted_gas_coords])

        dock_gas_coords = np.array([coord for coords in dock_gas_coords.values() for coord in coords])

        if adjusted_gas_coords.size == 0 or dock_gas_coords.size == 0:
            raise ValueError("One or both coordinate sets are empty.")

        # Calculate centroids
        centroid1 = np.mean(adjusted_gas_coords, axis=0)
        centroid2 = np.mean(dock_gas_coords, axis=0)

        # Center the geometries around their centroids
        centered_geom1 = adjusted_gas_coords - centroid1
        centered_geom2 = dock_gas_coords - centroid2

        # Compute the covariance matrix
        H = np.dot(centered_geom1.T, centered_geom2)

        # Singular Value Decomposition (SVD)
        U, _, Vt = np.linalg.svd(H)

        # Compute the rotation matrix
        rotation_matrix = np.dot(Vt.T, U.T)

        return rotation_matrix

    def rotate_vibdisps_vector(self, rotation_matrix):
        """
        Applies the rotation matrix to the vibrational displacement vectors to align them with the docked molecule's position.

        Args:
            imaginary_vibrational_displacements (np.ndarray): The vibrational displacement vectors from the gas-phase TS.
            rotation_matrix (np.ndarray): The 3x3 rotation matrix.

        Returns:
            np.ndarray: The rotated vibrational displacement vectors.

        Raises:
            ValueError: If the vibrational displacement vectors do not have 3 components or if the rotation matrix is not a 3x3 matrix.
        """
        if rotation_matrix.shape != (3, 3):
            raise ValueError("Rotation matrix must be a 3x3 matrix.")
        sorted_gas_data = sorted(self.gas_data, key=lambda x: x["element"])
        rotated_gas_data = []
        for entry in sorted_gas_data:
            transformed_vibdisps = np.dot(entry["vibdisps"], rotation_matrix.T)
            transformed_entry = entry.copy()  # Make a copy of the current entry
            transformed_entry["vibdisps"] = transformed_vibdisps  # Update 'vibdisps' with the transformed data
            rotated_gas_data.append(transformed_entry)

        return rotated_gas_data

    def create_modecar(self, args, molecular_oxygens, catalyst_hydrogens, gas_data):
        """
        Creates the MODECAR file with vibrational displacements for a dimer method calculation.

        This function adjusts the gas-phase coordinates to the unit cell, calculates the docking rotation matrix,
        and applies vibrational displacements to the solid structure based on the given molecular oxygens and
        catalyst hydrogens. The resulting vibrational modes are used to generate the MODECAR file, which
        defines the vibrational displacements for the dimer transition state search.

        Parameters:
        -----------
        args : argparse.Namespace
            Command-line arguments or configuration containing necessary parameters for docking and other operations.

        molecular_oxygens : list of int
            Indices of molecular oxygen atoms in the structure for inclusion in the MODECAR file.

        catalyst_hydrogens : list of int
            Indices of catalyst hydrogen atoms in the structure, whose modes will be excluded (set to small displacement).

        gas_data : list of dict
            A list of dictionaries containing gas-phase data for the molecular structure. Each dictionary contains
            information about atoms, including a flag `"O_docked"` to identify docked oxygen atoms.

        Returns:
        modes : list of list of float
            A list of vibrational displacements (x, y, z) for each atom in the structure. The modes for catalyst
            framework atoms are set to [0.001, 0.001, 0.001], while the modes for non-framework atoms are derived
            from the gas-phase vibrational displacements. Molecular oxygen atoms have their modes included, while
            catalyst hydrogens are excluded.

        Raises:
        IndexError:
            If the index of molecular oxygen or catalyst hydrogen indices exceeds the available atoms in the structure.
        """

        print("Adjusting gas-phase coordinates to unit cell...")
        gas_geom_and_vibdisps_to_unit_cell = self.adjust_to_unit_cell()
        print("Gas-phase coordinates adjusted.")
        print("Calculating docking rotation matrix...")
        docked_molecule_coords = self.get_docked_molecule_coords(args)
        print("Docking rotation matrix calculated.")
        docking_rotation = self.get_docking_rotations(gas_geom_and_vibdisps_to_unit_cell, docked_molecule_coords)

        modecar_vibdisps = self.rotate_vibdisps_vector(docking_rotation)

        # Prepare modes for MODECAR
        modes = []
        index, ox_index, h_index = 0, 0, 0
        o_docked_indexes = [i for i, entry in enumerate(gas_data) if entry["O_docked"]]

        for element, coords_list in self.solid_coords.items():
            for i in range(len(coords_list)):
                if element in CATALYST_FRAMEWORK_ELEMENTS:
                    modes.append([0.001, 0.001, 0.001])
                if element in CATALYST_NON_FRAMEWORK_ELEMENTS:
                    if not modecar_vibdisps[index]["O_docked"]:
                        modes.append(modecar_vibdisps[index]["vibdisps"])
                        index += 1

                if element == "O_docked":
                    if ox_index == len(molecular_oxygens):
                        break
                    for i in range(len(molecular_oxygens)):
                        modes.insert(
                            molecular_oxygens[i],
                            modecar_vibdisps[o_docked_indexes[i]]["vibdisps"],
                        )
                        print(
                            f"--include-molecular-oxygens, atom number {molecular_oxygens[i]} modes {modecar_vibdisps[o_docked_indexes[i]]['vibdisps']} written to output"
                        )
                        index += 1
                        ox_index += 1

                if element == "H_acid":
                    if h_index == len(catalyst_hydrogens):
                        break
                    for i in range(len(catalyst_hydrogens)):
                        modes[catalyst_hydrogens[i]] = [0.001, 0.001, 0.001]
                        print(
                            f"--exclude-catalyst-hydrogens atom number {catalyst_hydrogens[i]} modes set to [0.001, 0.001, 0.001] in output"
                        )
                        h_index += 1
                        index += 1

        return modes
