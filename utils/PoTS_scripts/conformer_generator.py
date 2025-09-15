#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Conformer generator script using RDKit.
Generates, minimizes, and clusters conformers from SMILES strings.
"""

import argparse
import os
import subprocess
import time
import random
from rdkit.Chem import MolFromSmiles, AddHs
from rdkit.Chem.AllChem import (
    EmbedMultipleConfs,
    UFFGetMoleculeForceField,
    MMFFGetMoleculeForceField,
    MMFFGetMoleculeProperties,
    GetConformerRMS,
)
from rdkit.Chem.rdmolops import RemoveHs

# Constants
UFF_ELEMENTS = ["B", "Al"]  # Switch to UFF if these elements are present
DEFAULT_GEOM_COMPARE_TIMEOUT = 300


# --------------------------
# Helper Functions
# --------------------------
def write_xyz(coords, filename, comment):
    """Write XYZ file from coordinates."""
    with open(filename, "w") as fp:
        fp.write(f"{len(coords)}\n")
        fp.write(f"{comment}\n")
        for atom in coords:
            fp.write(
                "%s %.4f %.4f %.4f\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2])
            )


def _extract_atomic_type(conformer):
    """Return list of atomic symbols for a conformer."""
    return [atom.GetSymbol() for atom in conformer.GetOwningMol().GetAtoms()]


def _atomic_pos_from_conformer(conformer):
    """Return list of atomic positions for a conformer."""
    return [
        [pos.x, pos.y, pos.z]
        for pos in (
            conformer.GetAtomPosition(i) for i in range(conformer.GetNumAtoms())
        )
    ]


# --------------------------
# Conformer Generator Class
# --------------------------
class ConformerGenerator:
    """Generates, minimizes, and clusters conformers for a molecule."""

    def __init__(self, smiles, forcefield="mmff"):
        self.mol = MolFromSmiles(smiles)
        self.full_clusters = []
        self.forcefield = forcefield
        self.conf_energies = []
        self.initial_confs = None
        self.smiles = smiles

    def generate(
        self,
        max_generated_conformers=50,
        prune_thresh=0.01,
        maxattempts_per_conformer=5,
        output=None,
        threads=1,
    ):
        """Generate conformers using RDKit."""
        self.mol = AddHs(self.mol, addCoords=True)
        self.initial_confs = EmbedMultipleConfs(
            self.mol,
            numConfs=max_generated_conformers,
            pruneRmsThresh=prune_thresh,
            maxAttempts=maxattempts_per_conformer,
            useRandomCoords=False,
            numThreads=threads,
            randomSeed=random.randint(1, 10000000),
        )

        if len(self.initial_confs) == 0:
            if output:
                output.write(
                    f"Retrying with {max_generated_conformers * 10} attempts and random coords\n"
                )
            self.initial_confs = EmbedMultipleConfs(
                self.mol,
                numConfs=max_generated_conformers,
                pruneRmsThresh=prune_thresh,
                useRandomCoords=True,
                maxAttempts=10 * maxattempts_per_conformer,
                numThreads=threads,
                randomSeed=random.randint(1, 10000000),
            )

        if output:
            output.write(f"Generated {len(self.initial_confs)} initial conformers\n")
        return self.initial_confs

    def minimise(self, output=None, minimize=True, fake_energies=False):
        """Minimize conformers using the specified force field."""
        if self.forcefield not in ["mmff", "uff"]:
            raise ValueError("Unrecognised force field")

        if self.forcefield == "mmff" and not fake_energies:
            props = MMFFGetMoleculeProperties(self.mol)
            for i in range(len(self.initial_confs)):
                potential = MMFFGetMoleculeForceField(self.mol, props, confId=i)
                if potential is None:
                    potential = UFFGetMoleculeForceField(self.mol, confId=i)
                if minimize:
                    potential.Minimize()
                self.conf_energies.append((i, potential.CalcEnergy()))
        else:
            for i in range(len(self.initial_confs)):
                if fake_energies:
                    energy = 1.5 * float(i)
                else:
                    potential = UFFGetMoleculeForceField(self.mol, confId=i)
                    if minimize:
                        potential.Minimize()
                    energy = potential.CalcEnergy()
                self.conf_energies.append((i, energy))

        self.conf_energies = sorted(self.conf_energies, key=lambda tup: tup[1])
        return self.mol

    def cluster(
        self,
        rms_tolerance=0.1,
        max_ranked_conformers=10,
        energy_window=5,
        Report_e_tol=10,
        output=None,
    ):
        """Cluster conformers based on RMSD and energy."""
        self.full_clusters = []
        self.mol_no_h = RemoveHs(self.mol)
        confs = self.conf_energies[:]
        ignore = []

        for i, (index_1, energy_1) in enumerate(confs):
            if i in ignore:
                continue
            clustered = [[self.mol.GetConformer(id=index_1), energy_1, 0.0]]
            ignore.append(i)
            for j, (index_2, energy_2) in enumerate(confs):
                if j in ignore:
                    continue
                if abs(energy_1 - energy_2) <= energy_window:
                    clustered.append([self.mol.GetConformer(id=index_2), energy_2, 0.0])
                    ignore.append(j)
            self.full_clusters.append(clustered)

        ranked_clusters = [
            cluster[0] for cluster in self.full_clusters[:max_ranked_conformers]
        ]
        return ranked_clusters


# --------------------------
# Main Workflow
# --------------------------
def main(
    smiles_file,
    forcefield="mmff",
    nconf=20,
    nconf_gen=200,
    maxattempts_per_conformer=5,
    e_window=5.0,
    rms_tol=0.1,
    prun_tol=0.01,
    log="confgen.log",
    rep_e_window=5.0,
    fallback_to_align=False,
    threads=1,
    output_dir=None,
):
    start_time = time.time()
    output_dir = output_dir or os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    with open(os.path.join(output_dir, log), "w") as output:
        # Load SMILES
        smiles = []
        with open(smiles_file, "r") as f:
            for line in f:
                smi = line.strip().split()[0]
                if "." in smi:
                    raise ValueError("Only single-molecule SMILES are supported")
                smiles.append(smi)
        output.write("SMILES strings to run:\n{}\n".format("\n".join(smiles)))

        for molecule in smiles:
            if any([el in molecule for el in UFF_ELEMENTS]):
                output.write("Switching to UFF (boron/aluminum present)\n")
                forcefield = "uff"

            output.write(f"Analysing SMILES: {molecule}\n")
            MolFromSmiles(molecule)

            # Generate conformers
            output.write("Generating initial conformers\n")
            confgen = ConformerGenerator(smiles=molecule, forcefield=forcefield)
            confgen.generate(
                max_generated_conformers=int(nconf_gen),
                maxattempts_per_conformer=int(maxattempts_per_conformer),
                prune_thresh=float(prun_tol),
                output=output,
                threads=threads,
            )

            # Minimise
            try:
                confgen.minimise(output=output)
            except RuntimeError as e:
                output.write(f"Minimization failed: {e}. Using fake energies.\n")
                confgen.minimise(output=output, fake_energies=True)

            output.write("Conformations minimised. Energies:\n")
            output.write("\n".join([str(en[1]) for en in confgen.conf_energies]) + "\n")

            # Cluster
            clustered_confs = confgen.cluster(
                rms_tolerance=float(rms_tol),
                max_ranked_conformers=int(nconf),
                energy_window=float(e_window),
                Report_e_tol=float(rep_e_window),
                output=output,
            )

            for i, conformer in enumerate(clustered_confs):
                pos = _atomic_pos_from_conformer(conformer[0])
                elements = _extract_atomic_type(conformer[0])
                coords = list(zip(elements, pos))
                xyz_file = os.path.join(output_dir, f"conf_{i + 1}.xyz")
                write_xyz(coords=coords, filename=xyz_file, comment=conformer[1])
                output.write(f"Cluster {i} energy: {conformer[1]}\n")


# --------------------------
# CLI
# --------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate and cluster conformers from SMILES"
    )
    parser.add_argument("smiles_file", help="File containing SMILES strings")
    parser.add_argument("--output_dir", help="Directory for output files")

    # Conformer generation parameters
    parser.add_argument(
        "-g",
        "--nconf_gen",
        type=int,
        default=200,
        help="Number of conformers to generate",
    )
    parser.add_argument(
        "-e", "--e_window", type=float, default=5.0, help="Energy window for clustering"
    )
    parser.add_argument(
        "-p", "--prun_tol", type=float, default=0.01, help="Prune RMS threshold"
    )
    parser.add_argument(
        "-E", "--rep_e_window", type=float, default=5.0, help="Report energy window"
    )
    parser.add_argument(
        "-t", "--rms_tol", type=float, default=0.1, help="RMS tolerance for clustering"
    )
    parser.add_argument(
        "-n",
        "--nconf",
        type=int,
        default=20,
        help="Number of conformers to keep after clustering",
    )
    parser.add_argument(
        "-m",
        "--maxattempts_per_conformer",
        type=int,
        default=5,
        help="Max embedding attempts per conformer",
    )
    parser.add_argument(
        "--fallback_to_align",
        action="store_true",
        help="Fallback to obabel --align if RMS fails",
    )

    args = parser.parse_args()

    main(
        smiles_file=args.smiles_file,
        output_dir=args.output_dir,
        nconf=args.nconf,
        nconf_gen=args.nconf_gen,
        maxattempts_per_conformer=args.maxattempts_per_conformer,
        e_window=args.e_window,
        rms_tol=args.rms_tol,
        prun_tol=args.prun_tol,
        rep_e_window=args.rep_e_window,
        fallback_to_align=args.fallback_to_align,
    )
