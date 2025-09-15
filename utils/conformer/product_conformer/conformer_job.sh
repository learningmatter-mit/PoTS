#!/bin/bash
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH -p sched_mit_rafagb
#SBATCH --mem-per-cpu=3800

cwd=$(pwd)

python conformer_generator.py \
    -g 1500 -e 0.30 -p 0.5 -E 3.5 -t 0.30 \
    -n 50 \
    -m 5 \
    --fallback_to_align smiles.smi

