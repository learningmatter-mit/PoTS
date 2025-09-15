PoTS: Pipeline for Transition State Search in Zeolite Catalysis

This repository accompanies the manuscript
DIMER-assisted automated transition state identification in zeolite catalysis (Nature Computational Science, 2025).

It provides a reproducible pipeline (PoTS.py) and supporting scripts for generating, optimizing, docking, and refining reaction intermediates and transition states across different simulation levels (XTB → ORCA → VASP).

PoTS/
├── PoTS.py              # Main pipeline script
├── smarts.txt           # Example SMARTS mapping file (reactant >> product)
├── smiles.txt           # Example SMILES input file
├── utils/               # Core utilities for the pipeline
│   ├── conformer/       # Conformer generation
│   ├── DFT_gas_phase/   # ORCA calculations (gas-phase opt, scan, TS)
│   ├── DFT_periodic_solid/ # VASP calculations (periodic optimizations, TS)
│   ├── docking/         # Zeolite docking and retrieval
│   └── retrieve_conformers.py
├── scripts/             # Internal modules for DIMER optimization
├── examples/            # Ready-to-run examples (POSCAR, CIF, modes)
└── template_jobs/       # Example job scripts (bash)


Quickstart

Prepare input molecules

Write SMILES strings in smiles.txt:

reactant: CCC1=C[CH+]C=CC1C(C)c1ccccc1
product:  CCC1-C=C[CH+]-C=C1C(C)c1ccccc1


Write the mapped SMARTS in smarts.txt:

[#6:0]...>>[#6:0]...


Run the pipeline

python PoTS.py --smiles-file smiles.txt --smarts-file smarts.txt


What happens under the hood

Conformer generation → utils/conformer/

Semiempirical XTB optimizations → xtb/

Gas-phase DFT optimizations (ORCA) → utils/DFT_gas_phase/orca_opt/

Docking into zeolite frameworks → utils/docking/

Periodic DFT refinements (VASP) → utils/DFT_periodic_solid/

TS search via scan + eigenvector following → utils/DFT_gas_phase/orca_scan/ and orca_TS/


SMARTS-based reaction parsing

utils/DFT_gas_phase/read_smarts.py extracts the moving hydrogen and connected carbons from the mapped SMARTS.

This information is passed to prepare_orca_scan.py to automatically build scan inputs like:

%geom
  Scan
    B <C_index> <H_index> = <initial_dist>, <final_dist>, <n_steps>
  end
end


Notes

Cluster submission (Slurm) is left to the user — this repo prepares all job folders and input files but does not submit them automatically.

ORCA and VASP binaries must be installed and accessible on your system.

For details of the algorithm and benchmarking, refer to the manuscript.