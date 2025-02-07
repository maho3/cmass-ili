#!/bin/bash
#SBATCH --job-name=charm_mtng  # Job name
#SBATCH --array=181,512,723,965         # Job array range for lhid
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=8            # Number of tasks
#SBATCH --time=24:00:00         # Time limit
#SBATCH --partition=gpu        # Partition name
#SBATCH --gpus-per-node=1       # Number of GPUs per node
#SBATCH --account=phy240043-gpu   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task

# SLURM_ARRAY_TASK_ID=150
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
baseoffset=2000

module restore cmass
conda activate cmassrun
lhid=$((SLURM_ARRAY_TASK_ID + baseoffset))

# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run

nbody=mtnglike
sim=fastpm
multisnapshot=True
extras="" # nbody.zf=0.500015"
L=3000
N=384
keys_to_check=(0.586220 0.606330 0.626440 0.646550 0.666660 0.686770 0.706880 0.726990 0.747100 0.767210)
# keys_to_check=()

outdir=/anvil/scratch/x-mho1/cmass-ili/$nbody/$sim/L$L-N$N
echo "outdir=$outdir"

# Loop through offsets and process files
for offset in 0; do  # for offset in $(seq 0 200 800); do
    loff=$((lhid + offset))
    postfix="nbody=$nbody sim=$sim nbody.lhid=$loff multisnapshot=$multisnapshot $extras"
    file=$outdir/$loff/halos.h5

    # Check if file exists
    if [ -f $file ]; then
        echo "File $file exists."
        all_keys_exist=true

        # Check if all keys exist in the file
        for key in $keys_to_check; do
            if ! h5ls $file | grep -q $key; then
                all_keys_exist=false
                echo "Key $key does not exist in $file."
                break
            fi
        done

        # Process file based on key existence
        if $all_keys_exist; then
            echo "All keys exist in $file. Skipping..."
        else
            echo "Not all keys exist in $file. Rerunning..."
            python -m cmass.bias.rho_to_halo $postfix
        fi
    else
        echo "File $file does not exist. Running..."
        python -m cmass.bias.rho_to_halo $postfix
    fi
done
