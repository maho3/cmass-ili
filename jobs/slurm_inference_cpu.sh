#!/bin/bash
#SBATCH --job-name=inference  # Job name
# #SBATCH --array=0-99  # Array range
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=4            # Number of tasks
#SBATCH --time=4:00:00         # Time limit
#SBATCH --partition=shared  # Partition name
#SBATCH --account=phy240043  # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task

# SLURM_ARRAY_TASK_ID=3
if [[ -z "$PS1" ]]; then
    export TQDM_DISABLE=0
fi

module restore cmass
conda activate cmassrun

# exp_index=null
net_index=$SLURM_ARRAY_TASK_ID

# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run

nbody=mtnglike
sim=fastpm
infer=default

halo=False
galaxy=False
ngc=True
sgc=False
mtng=False

extras="nbody.zf=0.500015"
device=cpu

suffix="nbody=$nbody sim=$sim infer=$infer infer.exp_index=$exp_index infer.net_index=$net_index"
suffix="$suffix infer.halo=$halo infer.galaxy=$galaxy"
suffix="$suffix infer.ngc_lightcone=$ngc infer.sgc_lightcone=$sgc infer.mtng_lightcone=$mtng"
suffix="$suffix infer.device=$device $extras"

echo "Running inference with $suffix"
# python -m cmass.infer.preprocess $suffix
# python -m cmass.infer.train $suffix net=tuning
python -m cmass.infer.validate $suffix
