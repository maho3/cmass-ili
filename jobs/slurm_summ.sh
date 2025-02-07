#!/bin/bash
#SBATCH --job-name=summ   # Job name
#SBATCH --array=0-999         # Job array range for lhid
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=8            # Number of tasks
#SBATCH --time=00:20:00         # Time limit
#SBATCH --partition=shared      # Partition name
#SBATCH --account=phy240043   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task

module restore cmass
source /anvil/scratch/x-mho1/anaconda3/bin/activate
conda activate cmassrun
# lhid=$SLURM_ARRAY_TASK_ID
# SLURM_ARRAY_TASK_ID=3


# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run


for offset in 0 1000; do
    lhid=$(($SLURM_ARRAY_TASK_ID+$offset))

    postfix="nbody=abacuslike sim=pinocchio nbody.N=1024 nbody.lhid=$lhid"

    python -m cmass.diagnostics.summ $postfix
done