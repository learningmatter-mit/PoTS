#!/usr/bin/env python3
"""
PoTS.py - Pipeline for automated transition state search in zeolite catalysis.

This script performs the following workflow:

1. Reads reactant and product SMILES from a file.
2. Optionally canonicalizes them using RDKit.
3. Prepares conformer jobs for reactant and product.
4. Launches conformer generation scripts and retrieves conformers.
5. Runs XTB optimization on conformers.
6. Prepares ORCA input folders for single-point or geometry optimizations.
7. Prepares docking folders for reactants, products, and TS structures.
8. Prepares VASP folders for periodic calculations.
9. Prepares ORCA scan input for reactive H detection from SMARTS.
10. Prepares ORCA eigenvector-following (TS) calculations.
11. Parses ORCA TS outputs and sets up VASP TS jobs.
12. Finally, runs the PoTS code to prepare a MODECAR file for a Dimer run with VASP.

Usage examples:
    python PoTS.py --smiles-file smiles.txt --smarts-file smarts.txt
"""


import os
import sys
import re
import shutil
import argparse
import subprocess


def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: 
            .smiles_file (str): Path to file containing reactant and product SMILES.
            .no_rdkit (bool): Skip RDKit validation if True.
            .smarts_file (str): Path to file containing reactant>>product SMARTS.
    """

    parser = argparse.ArgumentParser(
        description="PoTS pipeline launcher - parse reactant and product SMILES"
    )
    parser.add_argument(
        "--smiles-file",
        default="smiles.txt",
        help="Path to file containing reactant and product SMILES, one per line: reactant\\nproduct"
    )

    parser.add_argument(
        "--no-rdkit",
        action="store_true",
        help="Skip RDKit validation/canonicalization (useful if RDKit not installed)"
    )

    parser.add_argument(
        "--smarts-file",
        default="smarts.txt",
        help="Path to file containing reactant>>product SMARTS (default: smarts.txt)"
    )
    return parser.parse_args()

def canonicalize_smiles(smiles):
    """
    Canonicalize a SMILES string using RDKit.
    
    Args:
        smiles (str): SMILES string to canonicalize.
    
    Returns:
        tuple: (canonical_smiles (str) or None, error_msg (str) or None)
               error_msg is 'rdkit-not-available' or 'invalid-smiles' on failure.
    """

    try:
        from rdkit import Chem
    except Exception:
        return None, "rdkit-not-available"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "invalid-smiles"
    can = Chem.MolToSmiles(mol, isomericSmiles=True)
    return can, None


def parse_smarts_to_graph(smarts):
    """
    Convert a SMARTS string into an adjacency graph.
    
    Args:
        smarts (str): SMARTS string with atom mapping numbers.
    
    Returns:
        tuple:
            nodes (dict): {mapnum: atom_type}, e.g., {0:'6', 1:'1', ...}
            edges (dict): {mapnum: set(neighbor_mapnums)}, e.g., {0:{1,2}, ...}
    
    Handles branches using parentheses.
    """

    nodes = {}
    edges = {}
    stack = []
    last_atom = None

    atom_pat = re.compile(r'\[#(\d+)(?:\+?\d*)?:(\d+)\]')

    i = 0
    while i < len(smarts):
        char = smarts[i]
        if char == '(':
            stack.append(last_atom)
            i += 1
        elif char == ')':
            last_atom = stack.pop()
            i += 1
        elif char in '-=#~':
            i += 1
        elif char == '[':
            m = atom_pat.match(smarts, i)
            if m:
                atom_type, mapnum = m.groups()
                mapnum = int(mapnum)
                nodes[mapnum] = atom_type
                edges.setdefault(mapnum, set())
                if last_atom is not None:
                    edges[last_atom].add(mapnum)
                    edges[mapnum].add(last_atom)
                last_atom = mapnum
                i = m.end()
            else:
                i += 1
        else:
            i += 1
    return nodes, edges

def find_first_moving_h(react_smarts, prod_smarts):
    """
    Identify the first moving hydrogen atom between reactant and product.
    
    Args:
        react_smarts (str): SMARTS string for reactant.
        prod_smarts (str): SMARTS string for product.
    
    Returns:
        tuple: (h_index (int), c_index (int)) of moving hydrogen and its original carbon.
               Returns (None, None) if no moving H found.
    """

    react_nodes, react_edges = parse_smarts_to_graph(react_smarts)
    prod_nodes, prod_edges = parse_smarts_to_graph(prod_smarts)
    
    for atom_map, atom_type in react_nodes.items():
        if atom_type != '1':  # only H
            continue
        react_neighbors = react_edges[atom_map]
        prod_neighbors = prod_edges[atom_map]
        if react_neighbors != prod_neighbors:
            # first C neighbor in reactant
            react_parent = next(n for n in react_neighbors if react_nodes[n] != '1')
            # first C neighbor in product (optional, can return if needed)
            # prod_parent = next(n for n in prod_neighbors if prod_nodes[n] != '1')
            return atom_map, react_parent

    return None, None  # no moving H found


def prepare_and_launch_conformer_jobs(base_dir, reactant, product):
    """
    Prepare conformer generation jobs for reactant and product.
    
    Creates folders, copies template scripts, writes smiles.smi,
    and launches the conformer generation jobs.
    
    Args:
        base_dir (str): Path to base PoTS directory.
        reactant (str): Reactant SMILES.
        product (str): Product SMILES.
    
    Returns:
        None
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

    # Step 0: Read reactant and product SMILES from input file
    reactant = None
    product = None
    with open(args.smiles_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("reactant:"):
                reactant = line.split(":", 1)[1].strip()
            elif line.startswith("product:"):
                product = line.split(":", 1)[1].strip()

    if not reactant or not product:
        raise ValueError(f"Could not find reactant or product in {args.smiles_file}")
    
    # Step 0b: Canonicalize SMILES using RDKit, if available
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

    if not os.path.isfile(args.smarts_file):
        print(f"SMARTS file {args.smarts_file} not found.", file=sys.stderr)
        sys.exit(1)

    # Step 1: Prepare and launch conformer generation jobs for reactant and product
    base_dir = os.path.dirname(os.path.abspath(__file__))
    prepare_and_launch_conformer_jobs(base_dir, reactant, product)

    # Step 2: Retrieve generated conformers with external script
    xtb_dir = os.path.join(base_dir, "xtb")
    retriever = os.path.join(xtb_dir, "retrieve_conformers.py")

    for mol in ["reactant", "product"]:
        log_file = os.path.join(base_dir, "conformer", f"{mol}_conformer", "confgen.log")
        print(f"Running retriever for {mol}: {log_file}")
        subprocess.run(
            ["python", retriever, log_file],
            check=True
        )

    # Step 3: Launch XTB optimization jobs for each conformer
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

    # Step 4: Prepare ORCA input folders for optimized conformers
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

    # Step 6: Prepare VASP input folders for reactant and product optimizations
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
        
            
        print(f"Running create_vasp_opt_folders.py for {label}...")
        subprocess.run(
            ["python", vasp_script, label, docking_folder, vasp_opt_folder, incar_src, kpoints],
            cwd=vasp_opt_folder,
            check=True
        )
        
        # NOTE: VASP calculations themselves should be submitted to the cluster later

    
    # Step 7: Prepare ORCA scan inputs based on moving H detected from SMARTS
    smarts_file = os.path.join(base_dir, args.smarts_file)
    with open(smarts_file) as f:
        content = f.read().strip()

    # Expecting a single line: reactant>>product
    if '>>' not in content:
        print("SMARTS file must contain reactant>>product on one line")
        sys.exit(2)

    react_smarts, prod_smarts = content.split('>>')

    h_index, c_index = find_first_moving_h(react_smarts, prod_smarts)
    print('Identified reactant indexes from SMARTS', h_index, c_index)

    label = "product"  # or "product", depending on which one you want
    final_distance = 1.2  # target distance in Angstroms
    scan_script = os.path.join(base_dir, "DFT_gas_phase", "prepare_orca_scan.py")
    orca_scan_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_scan")

    # call the script
    subprocess.run(
        ["python", scan_script, label, str(h_index), str(c_index), str(final_distance)],
        check=True
    )

    print(f"ORCA scan input prepared for {label} in {orca_scan_folder}")

    # Step 8: Prepare ORCA TS (transition state) calculation folders
 
    print("Preparing ORCA EVF (transition state) calculation...")

    orca_scan_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_scan")
    orca_ts_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_TS")

    parse_script = os.path.join(base_dir, "DFT_gas_phase", "parse_scan_hei.py")

    subprocess.run(
        ["python", parse_script, orca_scan_folder, orca_ts_folder],
        cwd=orca_ts_folder,
        check=True
    )


    # Step 9: Prepare docking folders for TS structure
    print("Preparing docking folder for transition state...")

    orca_ts_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_TS")
    ts_xyz = os.path.join(orca_ts_folder, "orca_ts_opt.xyz")

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
        #job_sh = os.path.join(docking_ts_folder, "docking_job.sh")
        #if os.path.isfile(job_sh):
        #    subprocess.run(["bash", job_sh], cwd=docking_ts_folder, check=True)

    # Step 10: Prepare VASP folders for dimer TS calculation

    print("Preparing VASP folders for Dimer TS...")

    # ORCA optimized structures
    orca_folder = os.path.join(base_dir, "DFT_gas_phase", "orca_TS")
    # Corresponding docking folder
    docking_folder = os.path.join(base_dir, "docking", "ts_docking")
    # Destination VASP folder
    vasp_ts_folder = os.path.join(base_dir, "DFT_periodic_solid", "ts_vasp_dimer")
    os.makedirs(vasp_opt_folder, exist_ok=True)

    # Path to create_vasp_folders.py
    vasp_script = os.path.join(base_dir, "DFT_periodic_solid", "create_vasp_dimer_folder.py")
    if not os.path.isfile(vasp_script):
        print(f"Warning: {vasp_script} not found, skipping.")
        
    incar_src = os.path.join(base_dir, "DFT_periodic_solid", "INCAR-dimer")
    if not os.path.isfile(incar_src):
        print(f"Warning: {incar_src} not found, skipping.")
        
    kpoints = os.path.join(base_dir, "DFT_periodic_solid", "KPOINTS")
    if not os.path.isfile(kpoints):
        print(f"Warning: {kpoints} not found, skipping.")

    print(f"Running create_vasp_dimer_folders.py...")
    subprocess.run(
        ["python", vasp_script, docking_folder, vasp_ts_folder, incar_src, kpoints],
        cwd=vasp_opt_folder,
        check=True
    )

    print("Parsing ORCA TS output to vibdisps...")

    vasp_ts_folder = os.path.join(base_dir, "DFT_periodic_solid", "ts_vasp_dimer")
    orca_out = os.path.join(base_dir, "DFT_gas_phase", "orca_TS", "orca_ts.out")

    vib_script = os.path.join(base_dir, "DFT_periodic_solid", "parse_vibdisps.py")
    subprocess.run(
        ["python", vib_script, orca_out, vasp_ts_folder],
        check=True
    )

    # Step 11: Run PoTS code and prepare MODECAR for Dimer run with VASP

    # Path to the run script
    run_pots = os.path.expanduser("~/PoTS/scripts/run.py")

    target_dir = os.path.join(os.getcwd(), "DFT_periodic_solid", "ts_vasp_dimer")
    # Change directory
    os.chdir(target_dir)
    
    # Run the script
    subprocess.run(
        ["python", run_pots],
        check=True
    )

    # NOTE: VASP calculations themselves should be submitted to the cluster later


if __name__ == "__main__":
    main()

