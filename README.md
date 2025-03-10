# :link: DIMER: Transition State Mode Generation for Solid Catalysts

This repository provides scripts to generate vibrational displacement modes for DIMER transition state (TS) searches in crystalline structures.

This code has been tested after docking pose generation with [VOID](https://github.com/learningmatter-mit/VOID), but it should also work seamlessly with other pose generation software. 

## Overview

The code automates the generation of vibrational displacement modes for catalytic reactions, focusing on hydrocarbon reactions in acidic zeolites. The default output format is a `MODECAR` file for `VASP`, with atomic types in `POSCAR` and `POTCAR` files sorted alphabetically to ensure compatibility.

It identifies oxygen atoms in gas-phase molecules docked into the crystal and distinguishes them from those in the catalyst framework (e.g., zeolites or oxides). Hydrogens at catalyst acid sites are also recognized and excluded from reactant mode calculations. For more details and examples, see the [examples folder](https://github.com/learningmatter-mit/DIMER/tree/dimer_code/examples) or the end of this document.

The code is tailored for reactions in zeolites, primarily hydrocarbon conversions and some involving light heteroatoms. To handle a broader range of atom types or catalyst compositions, modify the `CATALYST_FRAMEWORK_ELEMENTS` and `CATALYST_NON_FRAMEWORK_ELEMENTS` variables in ([scripts/chemistry_composition.py](https://github.com/learningmatter-mit/DIMER/blob/dimer_code/scripts/chemistry_composition.py)), where additional reaction examples are also available.

Verbose output is enabled by default to aid users in understanding execution and debugging.

The code can also generate mode input files for other computational chemistry software, such as `CP2K`, `GULP`, `ASE`, and `LAMMPS`.

## Installation

To install this repository, follow these steps:

    git clone https://github.com/learningmatter-mit/DIMER

    pip install numpy

    pip install pymatgen


## Usage

To run the mode generation script, navigate to the directory containing the structure files. 

Make sure the following input files are available:

- A crystal structure file, either `POSCAR` or `.cif`. If both are available, the script prioritizes parsing the `POSCAR`.

- Gas molecule coordinates in an `.xyz` file.

- A `vibdisps` file containing negative-frequency TS gas-phase vibrational displacements.

- **Optional (but recommended):** A `POTCAR` file for `VASP` runs, which will be sorted alphabetically alongside the `POSCAR`.

Then execute:

    python ../../scripts/run.py

#### Optional Command-line Arguments

--cp2k: Generates the `&MODE_VECS` section for the `CP2K` dimer input file `CP2KMODEVECS`

--gulp: Generates the vector section for the `GULP` dimer input file `GULPvector`

--ase: Generates the mode vector section for the `ASE` dimer input file `ASE_INPUT`

--lammps: Generates the `displace_atoms` section for the LAMMPS dimer input file `LAMMPS_INPUT`

#### Using Optional Command-line Arguments

To generate initial guess modes for a DIMER TS search in formats other than `VASP`, add the appropriate arguments to the script. The script will write the required files to the folder, allowing you to request one or multiple formats simultaneously.

    python run.py --cp2k --gulp --ase --lammps
    
If no argument is provided, the script defaults to VASP format, with atomic types in `POSCAR` and `POTCAR` files sorted alphabetically.

## Examples

This section provides demonstrations of script usage across various scenarios and file types.

#### POSCAR_test (maybe turn this into LINKs when naming is clear)

This example uses a crystal structure in `POSCAR` format. The script automatically sorts atomic types alphabetically in `POSCAR` and `POTCAR` files, calculates the vibrational displacements, and generates a `MODECAR` file for `VASP`.

#### CIF_test

Analogous example to the `POSCAR_test` but using a `.cif` file for the crystal structure. 

#### Molecular_Oxygens_test

In this example, the docked gas-phase molecule contains oxygen atoms, which the script autonomously identifies and distinguishes from those in the catalyst framework. Vibrational displacements are then generated specifically for the gas-phase molecule.

#### Acid_Hydrogens_test

This example demonstrates the automatic detection and exclusion of hydrogen atoms from acidic catalyst sites during vibrational mode calculations, particularly for zeolite systems.

#### Molecular_Oxygens_Acid_Hydrogens_test

A combined example featuring a gas-phase molecule with oxygen atoms and a catalyst with acidic hydrogen atoms. The script detects and processes both the docked gas-phase oxygen atoms and the acid-site hydrogen atoms during vibrational mode calculations.

#### CP2K_ASE_GULP_LAMMPS_test

This example demonstrates how to request output formats for DIMER TS searches in software other than `VASP`, such as `CP2K`, `ASE`, `GULP`, and `LAMMPS`. Use the appropriate arguments (e.g., `--cp2k`, `--gulp`, `--ase`, `--lammps`) to generate the required input files for these software packages.

#### Modes_for_Metals

An example for generating vibrational displacement modes in metal-containing systems. For systems with a broader range of atom types or different catalyst compositions, adjust the variables in the [chemistry_composition.py](https://github.com/learningmatter-mit/DIMER/blob/dimer_code/scripts/chemistry_composition.py) file to include the specific metal elements. Refer to the `example_run.sh` script in [examples/Modes_for_Metals](https://github.com/learningmatter-mit/DIMER/tree/dimer_code/examples/Modes_for_Metals_test) for a detailed demonstration of the procedure.

## Citing
The publication describing the algorithm and the software is the following:

### LINK

If you use this software, please cite the paper above.
