#!/bin/bash
# Run XTB optimizations on all folders inside xtb_runs

# XTB binary
XTB_BIN="/home/paufv/miniforge3/envs/htvs/bin/xtb"

# Number of threads per run
NPROCS=4
OPT_LEVEL="normal"       # optimization level: normal, loose, verytight
TYPE_FLAG="--gfn2"       # or "--gfnff"
SOLVENT_FLAG=""           # e.g., "--alpb water" if needed

# Catch termination signals
trap 'kill -9 $pid; exit 1' SIGTERM SIGINT

# Loop over all folders

# Files
XYZ_FILE="gfn_opt.xyz"
OUT_FILE="gfn_opt.out"
INPUT_FILE="xtb.inp"

# Read charge and multiplicity from xtb.inp if you want to automate
CHARGE=0
MULTIPLICITY=1

# Run XTB
$XTB_BIN $XYZ_FILE \
    --opt $OPT_LEVEL \
    --parallel $NPROCS \
    --chrg $CHARGE \
    --uhf $((MULTIPLICITY - 1)) \
    $SOLVENT_FLAG $TYPE_FLAG \
    --input $INPUT_FILE > $OUT_FILE & pid=$!

wait $pid
cd - > /dev/null

