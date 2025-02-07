#!/bin/bash
#SBATCH --job-name=alternativehod   # Job name
#SBATCH --array=130-181         # Job array range for lhid
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=32            # Number of tasks
#SBATCH --time=00:30:00         # Time limit
#SBATCH --partition=shared      # Partition name
#SBATCH --account=phy240043   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task


module restore cmass
conda activate cmassrun

# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run

lhid=$SLURM_ARRAY_TASK_ID
Nhod=5

for i in $(seq 0 $(($Nhod-1))); do
    hod_seed=$((lhid*10+i+1))
    postfix="nbody=abacus sim=nbody-leauthaud nbody.lhid=$lhid bias.hod.seed=$hod_seed bias.hod.model=leauthaud11 bias.hod.default_params=behroozi10"
    echo "Running $postfix"
    python -m cmass.bias.apply_hod $postfix
    python -m cmass.diagnostics.summ $postfix diag.galaxy=True
done
