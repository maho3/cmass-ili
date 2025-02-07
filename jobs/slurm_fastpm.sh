#!/bin/bash
#SBATCH --job-name=quijotelike  # Job name
#SBATCH --array=100-199         # Job array range for lhid
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=128            # Number of tasks
#SBATCH --time=00:45:00         # Time limit
#SBATCH --partition=wholenode # Partition name
#SBATCH --account=phy240043   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task

# SLURM_ARRAY_TASK_ID=199
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
baseoffset=2000

module restore cmass
source /anvil/scratch/x-mho1/anaconda3/bin/activate
conda activate cmassrun
lhid=$((SLURM_ARRAY_TASK_ID + baseoffset))

# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run

nbody=quijotelike
sim=fastpm
extras="nbody.matchIC=0 meta.cosmofile=./params/big_sobol_params.txt"
L=1000
N=128

outdir=/anvil/scratch/x-mho1/cmass-ili/$nbody/$sim/L$L-N$N
echo "outdir=$outdir"

for offset in $(seq 0 200 1800); do
    loff=$((lhid + offset))
    
    postfix="nbody=$nbody sim=$sim nbody.lhid=$loff $extras"

    file=$outdir/$loff/nbody.h5
    if [ -f $file ]; then
        echo "File $file exists."
    else
        echo "File $file does not exist."
        python -m cmass.nbody.fastpm $postfix
    fi
done
