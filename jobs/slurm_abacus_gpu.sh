#!/bin/bash
#SBATCH --job-name=abacuspinn_gpu   # Job name
#SBATCH --array=0-999         # Job array range for lhid
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=32            # Number of tasks
#SBATCH --gpus-per-node=1     # Number of GPUs
#SBATCH --time=01:00:00         # Time limit
#SBATCH --partition=gpu      # Partition name
#SBATCH --account=phy240043-gpu   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task

# SLURM_ARRAY_TASK_ID=130

module restore cmassgpu
conda activate cmassgpu
lhid=$SLURM_ARRAY_TASK_ID


# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run

Nhod=5
Naug=1

nbody=pinocchio_2gpch
sim=pinocchio
multisnapshot=True
diag_from_scratch=True
rm_galaxies=False
extras="" # "nbody.zf=0.500015"
L=2000
N=1024

outdir=/anvil/scratch/x-mho1/cmass-ili/abacuslike/$sim/L$L-N$N
# outdir=/anvil/scratch/x-mho1/cmass-ili/$nbody/$sim/L$L-N$N
echo "outdir=$outdir"


for offset in 0 1000; do
    lhid=$(($SLURM_ARRAY_TASK_ID+offset))

    postfix="nbody=$nbody sim=$sim nbody.lhid=$lhid multisnapshot=$multisnapshot diag.from_scratch=$diag_from_scratch $extras"

    # density
    # python -m cmass.diagnostics.summ diag.density=True $postfix 

    # halos
    file=$outdir/$lhid/halos.h5
    if [ -f $file ]; then
        echo "File $file exists."
    else
        echo "File $file does not exist."
        # python -m cmass.bias.rho_to_halo $postfix
    fi
    python -m cmass.diagnostics.summ diag.halo=True $postfix 

    # # galaxies
    # for i in $(seq 0 $(($Nhod-1))); do
    #     hod_seed=$((lhid*10+i+1))
    #     printf -v hod_str "%05d" $hod_seed
    #     file=$outdir/$lhid/galaxies/hod$hod_str.h5
    #     if [ -f $file ]; then
    #         echo "File $file exists."
    #     else
    #         echo "File $file does not exist."
    #         python -m cmass.bias.apply_hod $postfix bias.hod.seed=$hod_seed
    #     fi
    #     python -m cmass.diagnostics.summ $postfix diag.galaxy=True bias.hod.seed=$hod_seed

    #     # augments
    #     for aug_seed in $(seq 0 $(($Naug-1))); do
    #         printf -v aug_str "%05d" $aug_seed
    #         # lightcone
    #         file=$outdir/$lhid/sgc_lightcone/hod${hod_str}_aug${aug_str}.h5
    #         if [ -f $file ]; then
    #             echo "File $file exists."
    #         else
    #             echo "File $file does not exist."
    #             python -m cmass.survey.selection survey=cmass_sgc $postfix bias.hod.seed=$hod_seed survey.aug_seed=$aug_seed
    #         fi
    #         python -m cmass.diagnostics.summ diag.sgc=True bias.hod.seed=$hod_seed survey.aug_seed=$aug_seed $postfix 
    #     done

    #     # Trash collection
    #     if [ $rm_galaxies = True ]; then
    #         # galaxies
    #         echo "Removing galaxies for lhid=$lhid hod_seed=$hod_seed"
    #         rm $outdir/$lhid/galaxies/hod$hod_str.h5

    #         # lightcone
    #         echo "Removing lightcone for lhid=$lhid hod_seed=$hod_seed"
    #         rm $outdir/$lhid/sgc_lightcone/hod${hod_str}_aug*.h5
    #     fi
    # done
done