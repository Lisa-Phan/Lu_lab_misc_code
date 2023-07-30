#!/bin/bash

#SBATCH -J testrosetta             # Job name
#SBATCH -o /scratch/09069/dhp563/slurm_input_output/out      # Name of stdout output file (%j expands to job ID)
#SBATCH -e /scratch/09069/dhp563/slurm_input_output/err      # Name of stderr output file (%j expands to job ID)
#SBATCH -p normal            # Queue name
#SBATCH -N 1                 # Total number of nodes
#SBATCH -n 3                 # Total number of mpi tasks
#SBATCH -t 3:00:00           # Run time (hh:mm:ss)

# Load any required modules or environment variables

python3 -m pipenv shell
python3 testrosetta.py RKR_Mb.pdb


