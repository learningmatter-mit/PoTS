"""
This module defines the composition of catalysts and reactants used in dimer calculations.
Depending on the specific reaction being studied, the `CATALYST_FRAMEWORK_ELEMENTS` and `CATALYST_NON_FRAMEWORK_ELEMENTS`
can be modified to reflect the appropriate elements.

The predefined elements are based on hydrocarbon reactions in acidic zeolites, but other catalyst compositions
used in dimer calculations are also provided as examples.

Refer to the README.md file for more details.
"""

CATALYST_FRAMEWORK_ELEMENTS = ['Al', 'Si', 'O']
CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H", "N", "F", "Cl", "Br", "I", "P", "S", "B"]

# Copper catalysed reactions in acidic zeolites where the Cu is fixed in the framework
# CATALYST_FRAMEWORK_ELEMENTS = ["Si", "O", "Al", "Cu"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H", "N"]

# Copper catalysed reactions in acidic zeolites where the Cu comes from the reactant
# CATALYST_FRAMEWORK_ELEMENTS = ["Si", "O", "Al"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H", "N", "Cu"]

# Platinum on Alumina Catalyzed Reactions
# CATALYST_FRAMEWORK_ELEMENTS = ["Al", "O"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["Pt", "H", "C"]

# Iron-Catalyzed Fischer-Tropsch Synthesis
# CATALYST_FRAMEWORK_ELEMENTS = ["Fe", "O"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H"]

# Zinc-Zirconium Catalyzed Methanol Synthesis
# CATALYST_FRAMEWORK_ELEMENTS = ["Zn", "Zr", "O"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H"]

# Gold Nanoparticle Catalyzed Oxidation
# CATALYST_FRAMEWORK_ELEMENTS = ["Au", "O", "Si"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H", "O"]

# Nickel-Catalyzed Methane Reforming
# CATALYST_FRAMEWORK_ELEMENTS = ["Ni", "O"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H"]

# Cobalt-Catalyzed Ammonia Synthesis
# CATALYST_FRAMEWORK_ELEMENTS = ["Co", "O"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["N", "H"]

# Vanadium-Catalyzed Olefin Metathesis
# CATALYST_FRAMEWORK_ELEMENTS = ["V", "O"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["C", "H"]

# Ruthenium-Catalyzed Water Splitting
# CATALYST_FRAMEWORK_ELEMENTS = ["Ru", "O", "Ti"]
# CATALYST_NON_FRAMEWORK_ELEMENTS = ["H", "O"]
