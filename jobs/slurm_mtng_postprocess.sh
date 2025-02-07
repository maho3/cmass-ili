#!/bin/bash
#SBATCH --job-name=mtnglike_postprocess   # Job name
#SBATCH --time=00:40:00         # Time limit
#SBATCH --account=phy240043   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=128            # Number of tasks
#SBATCH --partition=wholenode      # Partition name

# SLURM_ARRAY_TASK_ID=1457
# offset=2000
echo "Offset: $offset, SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

module restore cmass
conda activate cmassrun
lhid=$((SLURM_ARRAY_TASK_ID + offset))

cd /home/x-mho1/git/ltu-cmass-run
outdir=/anvil/scratch/x-mho1/cmass-ili/mtnglike/fastpm/L3000-N384

# check if nbody.h5 exists
file="$outdir/$lhid/nbody.h5"
if [ -f "$file" ]; then
    echo "File $file exists."
else
    echo "File $file does not exist."
    python -m cmass.nbody.postprocess_fastpm nbody=mtnglike nbody.lhid=$lhid
fi
