# PoTS — Pipeline for Transition-State Search 

PoTS (Python for Transition States) is a lightweight pipeline that automates the *local* workflow from molecular SMILES/SMARTS → conformers → semiempirical optimisation (xTB) → gas-phase DFT (ORCA) scan/TS preparation → docking → periodic VASP inputs.

Note: this pipeline driver (`PoTS.py`) runs from the `utils/` directory and works exclusively with the content in `utils/`.

---
## Software requirements

PoTS relies on Python and several external chemistry packages for different workflow steps. These must be installed separately:

- **RDKit** — required for conformer generation (`utils/conformer/conformer_generator.py`). Supports MMFF or UFF force fields.

- **Open Babel** (optional) — fallback for RMS alignment if clustering fails.

- **xTB** — semiempirical optimization (`xtb_job.sh`).

- **ORCA** — gas-phase DFT optimizations, scans, and TS calculations.

- [**VOID**](https://github.com/learningmatter-mit/VOID/tree/master?tab=readme-ov-file#void-voronoi-organic-inorganic-docker) — docking of molecules inside periodic frameworks accounting for the [cation-acid site](https://github.com/learningmatter-mit/VOID/tree/master?tab=readme-ov-file#docking-for-reactivity) distances (`utils/docking/`).

- **VASP** — periodic DFT calculations (geometry optimization and dimer TS search).

- **Slurm or local job submission system** — PoTS prepares job folders but does not execute DFT or xTB jobs.



## Quick overview

Run the pipeline (from the `utils/` directory):

```bash
cd /path/to/PoTS/utils
python PoTS.py --smiles-file smiles.txt --smarts-file smarts.txt
```

Or run from the repo root:

```
python utils/PoTS.py --smiles-file utils/smiles.txt --smarts-file utils/smarts.txt
```

- **smiles.txt** — file containing reactant / product SMILES (see *Input files*).

- **smarts.txt** — mapped SMARTS string (`reactant>>product`) used to detect the reactive indeces.

Important: PoTS.py automatically creates all the required folders and input files starting only from a smiles.txt and smarts.txt.
The utils/ folder in this repository already contains solved calculations and optimised files so new users can test the workflow from start to end without running heavy calculations. The adaptation of this pipeline to the user work environment is out of the reach of this tutorial folder.

---

### Input files (examples)

- **utils/smiles.txt** (plain text; two labelled lines — `reactant:` and `product:`):

**reactant**: CCC1=C[CH+]C=CC1C(C)c1ccccc1

**product**:  CCC1-C=C[CH+]-C=C1C(C)c1ccccc1

- **utils/smarts.txt** (single-line mapped SMARTS `reactant>>product`):
[#6:0](-[#1:17])1=[#6:1]([#1:18])...>>[#6:0](-[#1:17])1-[#6+1:1]...

(SMARTS must include atom-map indices `:N` so PoTS can detect reactant H and C atoms.)

---

**What PoTS does (high-level pipeline)**

PoTS orchestrates the following steps, operating on utils/ contents and creating structured folders under utils/:

- `Read inputs` — read smiles.txt and smarts.txt (both in utils/).

- (Optional) Canonicalize SMILES — use RDKit if available.

- `Conformer generation` — prepares utils/conformer/reactant_conformer and utils/conformer/product_conformer, copies conformer_generator.py and conformer_job.sh, writes smiles.smi, and launches conformer jobs.

- `Retrieve conformers` — runs utils/retrieve_conformers.py (or your configured retrieval script) to prepare XTB run folders.

- `XTB optimizations` — launches xtb_job.sh inside utils/xtb/<reactant|product>_xtb/conf_* folders.

- `Prepare ORCA inputs for minima` — utils/DFT_gas_phase/orca_opt/* are created from the best XTB outputs (script: utils/DFT_gas_phase/retrieve_xtbs.py).

- `Docking` — utils/docking/<reactant|product>_docking are prepared (script: utils/docking/retrieve_mins_to_dock.py). Docking jobs are not automatically submitted by PoTS by default.

- `Create VASP input folders` — utils/DFT_periodic_solid/<reactant|product>_vasp_opt created using create_vasp_opt_folders.py.

- `Prepare ORCA scan for TS search` — use mapped SMARTS to identify the moving hydrogen (script in utils/), compute initial distance from the chosen ORCA geometry, and write %geom Scan input (utils/DFT_gas_phase/prepare_orca_scan.py).

- `Run eigenvector-following / TS setup` — parse scan outputs, prepare EVF inputs and ORCA TS folders (utils/DFT_gas_phase/parse_scan_hei.py → orca_TS).

- `TS docking and VASP dimer` — prepare utils/docking/ts_docking and utils/DFT_periodic_solid/ts_vasp_dimer.

---

**Where to look (key files / locations)**

All paths below are relative to the repository root.

- `utils/PoTS.py` — main pipeline driver (run from utils/).

- `utils/smiles.txt` — SMILES input (reactant/product).

- `utils/smarts.txt` — mapped SMARTS (single-line reactant>>product).

- `utils/conformer/` — conformer generator scripts & template job scripts:
  - `conformer_generator.py`, `conformer_job.sh`

- `utils/xtb/` — folders for XTB runs (`reactant_xtb`, `product_xtb`) and template `xtb_job.sh`.

- `utils/DFT_gas_phase/` — ORCA inputs and helpers:
  - `orca_opt/` — ORCA input folders for minima.
  - `orca_scan/` — scan files, scan outputs.
  - `orca_TS/` — TS EVF inputs & outputs.
  - `prepare_orca_scan.py` — create scan input from ORCA geometry and SMARTS indices.
  - `parse_scan_hei.py` — parse scan outputs to EVF inputs.

- `utils/docking/` — docking templates and retrieval script `retrieve_mins_to_dock.py`.

- `utils/DFT_periodic_solid/` — VASP preparers:
  - `create_vasp_opt_folders.py`
  - `create_vasp_dimer_folder.py`
  - `parse_vibdisps.py`
  - target folders: `reactant_vasp_opt`, `product_vasp_opt`, `ts_vasp_dimer`

- `scripts/run.py` — optional final PoTS runner (invoked from `ts_vasp_dimer`).

---

**During the run PoTS will:**

- create `utils/conformer/*` jobs and launch them
- generate `utils/xtb/*` folders
- prepare `utils/DFT_gas_phase/*` and `utils/docking/*` inputs
- prepare `utils/DFT_periodic_solid/*` VASP folders
- create `utils/DFT_gas_phase/orca_scan/orca_scan.inp` (scan configuration) and `orca_TS` inputs

