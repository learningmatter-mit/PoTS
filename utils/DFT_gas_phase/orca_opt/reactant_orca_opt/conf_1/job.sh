#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 4300
#SBATCH -p sched_mit_rafagb
#SBATCH --mem-per-cpu=1900
#SBATCH --no-requeue
#SBATCH --signal=B:2@300

source ~/.bashrc

# import set_local_scratch, make_scratch_folder, run_orca, and clean_up functions
source ${HTVSDIR}/chemconfigs/orca/defaults/utils/orca_bash_funcs.sh

# set variables used by default functions
PLATFORM=engaging
MPI=/nfs/rafagblab001/software/mpi/openmpi313 # MPI module to load on supercloud or path to MPI for ORCA on engaging
ORCA=/nfs/rafagblab001/software/orca/orca_4_1_1_linux_x86-64_openmpi313

cwd=$(pwd)
set_local_scratch
make_scratch_folder

shopt -s extglob
trap 'kill -9 $pid; clean_up;  exit 1' SIGTERM SIGINT

run_orca
wait
clean_up

# LINES BELOW ARE ADDED AUTOMATICALLY BY JOB MANAGER #
touch job_manager-complete

