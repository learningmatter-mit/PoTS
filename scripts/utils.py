import os
from math import cos, sin, sqrt, radians
from chemistry_composition import CATALYST_FRAMEWORK_ELEMENTS


class DimerUtils:
    def convert_lattice_parameters_to_vectors(a, b, c, alpha, beta, gamma):
        """
        Convert lattice parameters to lattice vectors.

        Args:
            a (float): Lattice parameter a.
            b (float): Lattice parameter b.
            c (float): Lattice parameter c.
            alpha (float): Angle between b and c (in degrees).
            beta (float): Angle between a and c (in degrees).
            gamma (float): Angle between a and b (in degrees).

        Returns:
            list: A list of lattice vectors.

        Raises:
            ValueError: If the lattice parameters or angles are invalid, causing a math domain error.
        """
        try:
            # Convert angles from degrees to radians
            alpha_rad = radians(alpha)
            beta_rad = radians(beta)
            gamma_rad = radians(gamma)

            # Calculate the components of the lattice vectors
            v_x = a
            v_y = b * cos(gamma_rad)
            v_z = c * cos(beta_rad)
            v_yz = b * sin(gamma_rad)

            # Calculate the value inside the square root
            inside_sqrt = c**2 - v_z**2 - v_yz**2

            # Check if the value inside the square root is non-negative
            if inside_sqrt < 0:
                raise ValueError(
                    f"Math domain error: Attempted to calculate sqrt of a negative value. "
                    f"Check lattice parameters or angles. Calculated inside_sqrt = {inside_sqrt}, "
                    f"with c = {c}, v_z = {v_z}, v_yz = {v_yz}."
                )

            # Compute the third component of the lattice vector
            v_z = sqrt(inside_sqrt)

            # Construct the lattice vectors
            lattice_vectors = [
                [v_x, 0, 0],
                [v_y, v_yz, 0],
                [
                    v_z,
                    v_yz * (cos(alpha_rad) - cos(beta_rad) * cos(gamma_rad)) / sin(gamma_rad),
                    v_z,
                ],
            ]

            return lattice_vectors

        except ValueError as e:
            raise ValueError(
                f"Error in convert_lattice_parameters_to_vectors: {str(e)}. "
                "This typically occurs if the provided lattice parameters or angles "
                "result in invalid geometric calculations."
            )
        except Exception as e:
            raise Exception(
                f"Unexpected error in convert_lattice_parameters_to_vectors: {str(e)}. "
                "Please check the input parameters for any anomalies."
            )

    def direct_to_cartesian(coord, lattice_vectors):
        """
        Converts direct (fractional) coordinates to Cartesian coordinates.

        Args:
            coord (list): Direct coordinates [x, y, z].
            lattice_vectors (list): 3x3 list representing the lattice vectors.

        Returns:
            list: Cartesian coordinates [x, y, z].
        """
        cartesian_coord = [0.0, 0.0, 0.0]

        for i in range(3):
            for j in range(3):
                cartesian_coord[i] += coord[j] * lattice_vectors[j][i]

        return cartesian_coord

    def get_molecular_oxygen_catalyst_hydrogens_indicies(sorted_structure, threshold_O=4, threshold_H=1.2, bond_threshold_framework=1.7):
        """
        Check if any oxygen atom ("O") is close enough to an element different than "Si", "Al", or "O".
        Works with a pymatgen Structure object.

        Args:
            sorted_structure (Structure): A pymatgen Structure object.
            threshold_O (float): Distance threshold for oxygen atoms to check for proximity.
            threshold_H (float): Distance threshold for hydrogen atoms to check for proximity.
            bond_threshold_framework (float): Distance cutoff to exclude O near framework atoms (e.g., Si, Al).

        Returns:
            list: Sorted list of oxygen indices.
            list: Sorted list of hydrogen indices.
        """

        oxygen_indexes, hydrogen_indexes = [], []

        # Iterate over all pairs of atoms in the structure and calculate distances
        for i, site1 in enumerate(sorted_structure):
            element1 = site1.species_string  # Get the species symbol of site1

            for j, site2 in enumerate(sorted_structure):
                element2 = site2.species_string  # Get the species symbol of site2
                distance = sorted_structure.get_distance(i, j)  # Get the distance between site1 and site2

                # Check for hydrogen and catalyst framework elements based on the threshold_H
                if element1 == "H" and element2 in CATALYST_FRAMEWORK_ELEMENTS and j not in oxygen_indexes and distance < threshold_H:
                    if i not in hydrogen_indexes:
                        hydrogen_indexes.append(i)
                elif element2 == "H" and element1 in CATALYST_FRAMEWORK_ELEMENTS and i not in oxygen_indexes and distance < threshold_H:
                    if j not in hydrogen_indexes:
                        hydrogen_indexes.append(j)

                # Check for oxygen proximity based on the threshold_O
                if element1 == "O":
                    if element2 not in CATALYST_FRAMEWORK_ELEMENTS and distance < threshold_O:
                    # Check if O is also near a framework atom
                        near_framework = any(
                            sorted_structure.get_distance(i, k) < bond_threshold_framework and
                            sorted_structure[k].species_string in CATALYST_FRAMEWORK_ELEMENTS
                            for k in range(len(sorted_structure)) if k != i
                        )
                        if not near_framework and i not in hydrogen_indexes and i not in oxygen_indexes:
                            oxygen_indexes.append(i)
                elif element2 == "O":
                    if element1 not in CATALYST_FRAMEWORK_ELEMENTS and distance < threshold_O:
                        near_framework = any(
                            sorted_structure.get_distance(j, k) < bond_threshold_framework and
                            sorted_structure[k].species_string in CATALYST_FRAMEWORK_ELEMENTS
                            for k in range(len(sorted_structure)) if k != j
                        )
                        if not near_framework and j not in hydrogen_indexes and j not in oxygen_indexes:
                            oxygen_indexes.append(j)


        return sorted(oxygen_indexes), sorted(hydrogen_indexes)

    def indentify_molecular_oxygens_catalyst_hydrogens(structure):
        """
        Identifies molecular oxygen and catalyst hydrogen atoms in a given structure.

        This function sorts the input structure based on the element species and identifies
        the indices of molecular oxygen (O) and catalyst hydrogen (H) atoms. It performs validation
        to ensure that the identified atoms are within the bounds of the structure and match the
        species provided in the input structure.

        Parameters:
        structure : pymatgen.core.structure.Structure
            A pymatgen Structure object containing atomic coordinates and species data.

        Returns:
        molecular_oxygens : list of int
            A list of indices corresponding to molecular oxygen atoms in the structure.

        catalyst_hydrogens : list of int
            A list of indices corresponding to catalyst hydrogen atoms in the structure.

        Raises:
        ValueError:
            If any of the provided oxygen or hydrogen indices are out of bounds or do not match the
            expected species in the structure.
        """
        sorted_structure = structure.get_sorted_structure(key=lambda s: s.species_string)

        molecular_oxygens, catalyst_hydrogens = DimerUtils.get_molecular_oxygen_catalyst_hydrogens_indicies(sorted_structure)

        # Get the total number of atoms in the structure
        num_atoms = len(sorted_structure)

        # Initialize variables
        current_index = 0
        oxygen_indices, hydrogen_indices = [], []

        # Iterate over all atoms in the pymatgen structure
        for site in sorted_structure:
            element = site.species_string  # Get the element symbol (e.g., "O", "H", etc.)

            # Check if the element is "O" or "H" and collect their indices
            if element == "O":
                oxygen_indices.append(current_index)
            elif element == "H":
                hydrogen_indices.append(current_index)

            # Increment the current index to keep track of global positions
            current_index += 1

        if any(i < 1 or i > num_atoms for i in molecular_oxygens):
            raise ValueError("One or more provided oxygen indices are out of bounds.")
        if any(i not in oxygen_indices for i in molecular_oxygens):
            raise ValueError("One or more provided oxygen indices do not belong to indexed O atoms from the input file.")

        if any(i < 1 or i > num_atoms for i in catalyst_hydrogens):
            raise ValueError("One or more provided hydrogens indices are out of bounds.")
        if any(i not in hydrogen_indices for i in catalyst_hydrogens):
            raise ValueError("One or more provided hydrogen indices do not belong to indexed H atoms from the input file.")
        print(f"Include molecular oxygens: {molecular_oxygens}")
        print(f"Exclude catalyst hydrogens: {catalyst_hydrogens}")

        return molecular_oxygens, catalyst_hydrogens

    def is_file_empty(filepath):
        return os.stat(filepath).st_size == 0
