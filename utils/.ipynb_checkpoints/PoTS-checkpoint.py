#!/usr/bin/env python3
"""
PoTS.py

Usage examples:
    python PoTS.py "CCO" "CC=O"
    python PoTS.py --no-rdkit "CCO" "CC=O"
"""

import os
import sys
import shutil
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(
        description="PoTS pipeline launcher - parse reactant and product SMILES"
    )
    parser.add_argument(
        "reactant_smiles",
        help="Reactant SMILES string (wrap in quotes if it contains special chars)"
    )
    parser.add_argument(
        "product_smiles",
        help="Product SMILES string (wrap in quotes if it contains special chars)"
    )
    parser.add_argument(
        "--no-rdkit",
        action="store_true",
        help="Skip RDKit validation/canonicalization (useful if RDKit not installed)"
    )
    return parser.parse_args()

def canonicalize_smiles(smiles):
    """Return canonical SMILES using RDKit or (None, error_msg) on failure."""
    try:
        from rdkit import Chem
    except Exception:
        return None, "rdkit-not-available"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "invalid-smiles"
    can = Chem.MolToSmiles(mol, isomericSmiles=True)
    return can, None

def prepare_and_launch_conformer_jobs(base_dir, reactant, product):
    """
    Prepare reactant_ and product_ conformer job folders,
    copy job scripts, write smiles.smi, and launch the jobs.
    """
    conformer_dir = os.path.join(base_dir, "conformer")
    template_files = [
        os.path.join(conformer_dir, "conformer_generator.py"),
        os.path.join(conformer_dir, "conformer_job.sh")
    ]

    jobs = {
        "reactant_conformer": reactant,
        "product_conformer": product
    }

    for job_name, smiles in jobs.items():
        job_path = os.path.join(conformer_dir, job_name)
        os.makedirs(job_path, exist_ok=True)

        # Copy template files
        for f in template_files:
            shutil.copy(f, job_path)

        # Write smiles.smi
        smiles_file = os.path.join(job_path, "smiles.smi")
        with open(smiles_file, "w") as f:
            f.write(smiles + "\n")

        print(f"Prepared {job_name} at {job_path}")

        # Launch conformer job
        subprocess.run(["bash", "conformer_job.sh"], cwd=job_path, check=True)
        print(f"Launched conformer job in {job_path}")

def main():
    args = parse_args()
    reactant = args.reactant_smiles
    product = args.product_smiles

    if not args.no_rdkit:
        for name, smi in (("reactant", reactant), ("product", product)):
            can, err = canonicalize_smiles(smi)
            if err == "rdkit-not-available":
                print("RDKit not available: skipping SMILES validation/canonicalization.")
                break
            if err == "invalid-smiles":
                print(f"ERROR: {name} SMILES is invalid: {smi}")
                sys.exit(2)
            if name == "reactant":
                reactant = can
            else:
                product = can

    print("Parsed inputs:")
    print("  Reactant SMILES:", reactant)
    print("  Product  SMILES:", product)

    # Step 1: prepare and launch conformer jobs
    base_dir = os.path.dirname(os.path.abspath(__file__))
    prepare_and_launch_conformer_jobs(base_dir, reactant, product)

    # Step 2: retrieve conformers with external script
    xtb_dir = os.path.join(base_dir, "xtb")
    retriever = os.path.join(xtb_dir, "retrieve_conformers.py")

    for mol in ["reactant", "product"]:
        log_file = os.path.join(base_dir, "conformer", f"{mol}_conformer", "confgen.log")
        print(f"Running retriever for {mol}: {log_file}")
        subprocess.run(
            ["python", retriever, log_file],
            check=True
        )

    # Step 3: run XTB jobs
    for label in ["reactant", "product"]:
        xtb_dir = os.path.join(base_dir,"xtb", f"{label}_xtb")
        if not os.path.isdir(xtb_dir):
            print(f"Warning: {xtb_dir} not found, skipping XTB step for {label}")
            continue

        for conf in sorted(os.listdir(xtb_dir)):
            conf_path = os.path.join(xtb_dir, conf)
            job_script = os.path.join(conf_path, "xtb_job.sh")
            if os.path.isfile(job_script):
                print(f"Launching XTB job: {job_script}")
                subprocess.run(["bash", job_script], cwd=conf_path, check=True)

    # Step 4: Prepare ORCA opts for mins
    print("Preparing ORCA inputs (not submitting jobs)...")

    for label in ["reactant", "product"]:
        retriever = os.path.join(base_dir, "DFT_gas_phase", "orca_opt", "retrieve_xtbs.py")
        subprocess.run(
            ["python", retriever, label],
            check=True
        )


    # NOTE for users: ORCA jobs must be launched manually on Slurm clusters

    # Step 5: Parse Orca Optimized intermediates structures and run a normal docking cation-anion with VOID
    print("Preparing docking folders for reactants and products...")
    
    for label in ["reactant", "product"]:
        orca_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_opt", f"{label}_orca_opt")
        docking_folder = os.path.join(base_dir, "docking", f"{label}_docking")
        os.makedirs(docking_folder, exist_ok=True)
    
        # Copy extra files from the original docking/ folder
        docking_source_folder = os.path.join(base_dir, "docking")
        for extra in ["zeolite.cif", "docking_job.sh"]:
            src_file = os.path.join(docking_source_folder, extra)
            dst_file = os.path.join(docking_folder, extra)
            if os.path.isfile(src_file):
                shutil.copy2(src_file, dst_file)
            else:
                print(f"Warning: {src_file} not found, skipping.")
    
        retriever = os.path.join(docking_source_folder, "retrieve_mins_to_dock.py")
        
        # Run the retriever script in the docking folder, passing label and ORCA path
        subprocess.run(
            ["python", retriever, label, orca_folder],
            cwd=docking_folder,
            check=True
        )

    
        # Optionally launch docking_job.sh in each conf_ folder
        #for conf in sorted(os.listdir(docking_folder)):
        #    conf_path = os.path.join(docking_folder, conf)
        #    job_sh = os.path.join(conf_path, "docking_job.sh")
        #    if os.path.isdir(conf_path) and os.path.isfile(job_sh):
        #        subprocess.run(["bash", job_sh], cwd=conf_path, check=True)

    # Step 6: Prepare VASP folders for reactants and products
    print("Preparing VASP folders for reactants and products...")

    for label in ["reactant", "product"]:
        # ORCA optimized structures
        orca_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_opt", f"{label}_orca_opt")
        # Corresponding docking folder
        docking_folder = os.path.join(base_dir, "docking", f"{label}_docking")
        # Destination VASP folder
        vasp_opt_folder = os.path.join(base_dir, "DFT_periodic_solid", f"{label}_vasp_opt")
        os.makedirs(vasp_opt_folder, exist_ok=True)
    
        # Path to create_vasp_folders.py
        vasp_script = os.path.join(base_dir, "DFT_periodic_solid", "create_vasp_opt_folders.py")
        if not os.path.isfile(vasp_script):
            print(f"Warning: {vasp_script} not found, skipping.")
            continue

        incar_src = os.path.join(base_dir, "DFT_periodic_solid", "INCAR-opt")
        if not os.path.isfile(incar_src):
            print(f"Warning: {incar_src} not found, skipping.")
            continue

        kpoints = os.path.join(base_dir, "DFT_periodic_solid", "KPOINTS")
        if not os.path.isfile(kpoints):
            print(f"Warning: {kpoints} not found, skipping.")
            continue
        
            
        print(f"Running create_vasp_folders.py for {label}...")
        subprocess.run(
            ["python", vasp_script, label, docking_folder, vasp_opt_folder, incar_src, kpoints],
            cwd=vasp_opt_folder,
            check=True
        )
        
        # NOTE: VASP calculations themselves should be submitted to the cluster later



    # Step 7, Reactivity, Orca Scan


    # Step 8, From Orca Scan to Orca Eigenvector Following

 
    print("Preparing ORCA EVF (transition state) calculation...")

    orca_scan_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_scan")
    orca_ts_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_TS")

    parse_script = os.path.join(base_dir, "DFT_gas_phase", "parse_scan_hei.py")

    subprocess.run(
        ["python", parse_script, orca_scan_folder, orca_ts_folder],
        cwd=orca_ts_folder,
        check=True
    )

    # Dock the evf final TS geometry
    ## ADD CHECKING FOR NEGFREQS

    # Step 9: Parse Orca TS structure and prepare docking
    print("Preparing docking folder for transition state...")

    orca_ts_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_TS")
    ts_xyz = os.path.join(orca_ts_folder, "orca_opt.xyz")

    if not os.path.isfile(ts_xyz):
        print(f"Warning: {ts_xyz} not found, skipping TS docking.")
    else:
        docking_ts_folder = os.path.join(base_dir, "docking", "ts_docking", "conf_1")
        os.makedirs(docking_ts_folder, exist_ok=True)

        # Copy TS geometry as molecule.xyz
        shutil.copy2(ts_xyz, os.path.join(docking_ts_folder, "molecule.xyz"))
        print(f"Copied TS geometry {ts_xyz} → {docking_ts_folder}/molecule.xyz")

        # Copy zeolite.cif and docking_job.sh
        docking_source_folder = os.path.join(base_dir, "docking")
        for extra in ["zeolite.cif", "docking_job.sh"]:
            src_file = os.path.join(docking_source_folder, extra)
            dst_file = os.path.join(docking_ts_folder, extra)
            if os.path.isfile(src_file):
                shutil.copy2(src_file, dst_file)
                print(f"Copied {src_file} → {dst_file}")
            else:
                print(f"Warning: {src_file} not found, skipping.")

        # Optionally launch TS docking
        job_sh = os.path.join(docking_ts_folder, "docking_job.sh")
        if os.path.isfile(job_sh):
            subprocess.run(["bash", job_sh], cwd=docking_ts_folder, check=True)

    # Step 10: PoTS

    # Step 11: Dimer




if __name__ == "__main__":
    main()

