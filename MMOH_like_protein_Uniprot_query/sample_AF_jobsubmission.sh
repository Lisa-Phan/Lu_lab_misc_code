#!/bin/bash
# -----------------------------------------------------------------
#SBATCH -J my_af2_job                 # Job name
#SBATCH -o my_af2_job.%j.out          # Name of stdout output file
#SBATCH -e my_af2_job.%j.err          # Name of stderr output file
#SBATCH -p rtx                        # Queue (partition) name
#SBATCH -N 1                          # Total # of nodes
#SBATCH -n 1                          # Total # of mpi tasks
#SBATCH -t 10:00:00                   # Run time (hh:mm:ss)
#SBATCH -A myproject                  # Project/Allocation name
# -----------------------------------------------------------------

# Load modules (example path on ls6)
module unload xalt
module use /scratch/tacc/apps/bio/alphafold/modulefiles
module load alphafold/2.2.0-ctr

# Run AlphaFold
run_alphafold.sh --flagfile=$AF2_HOME/test/flags/full_dbs.ff \
                 --fasta_paths=$SCRATCH/ESM_fold_setup_test/esm/examples/protein-programming-language/mini.fasta \
                 --output_dir=$SCRATCH/ESM_fold_setup_test/esm/examples/protein-programming-language/Alphafold_prediction_mini \
                 --model_preset=monomer \
                 --use_gpu_relax=True