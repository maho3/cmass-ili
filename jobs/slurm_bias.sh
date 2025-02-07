#!/bin/bash
#SBATCH --job-name=quijotelike_bias   # Job name
#SBATCH --array=0-999         # Job array range for lhid
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=16            # Number of tasks
#SBATCH --time=03:00:00         # Time limit
#SBATCH --partition=shared      # Partition name
#SBATCH --account=phy240043   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task

# SLURM_ARRAY_TASK_ID=3

module restore cmass
source /anvil/scratch/x-mho1/anaconda3/bin/activate
conda activate cmassrun
lhid=$SLURM_ARRAY_TASK_ID

outdir=/anvil/scratch/x-mho1/cmass-ili/quijotelike/fastpm/L1000-N128

# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run

Nhod=10
Naug=5


for offset in 0 1000; do
    lhid=$(($SLURM_ARRAY_TASK_ID+offset))

    postfix="nbody=quijotelike sim=fastpm nbody.lhid=$lhid"

    # halos
    file=$outdir/$lhid/halos.h5
    if [ -f $file ]; then
        echo "File $file exists."
    else
        echo "File $file does not exist."
        python -m cmass.bias.rho_to_halo $postfix
    fi
    # galaxies
    for i in $(seq 0 $Nhod); do
        hod_seed=$((lhid*10+i+1))
        printf -v hod_str "%05d" $hod_seed
        file=$outdir/$lhid/galaxies/hod$hod_str.h5
        if [ -f $file ]; then
            echo "File $file exists."
        else
            echo "File $file does not exist."
            python -m cmass.bias.apply_hod $postfix bias.hod.seed=$hod_seed
        fi

        # # augments
        # for aug_seed in $(seq 0 $Naug); do
        #     printf -v aug_str "%05d" $aug_seed
        #     # lightcone
        #     file=$outdir/$lhid/lightcone/hod${hod_str}_aug${aug_str}.h5
        #     if [ -f $file ]; then
        #         echo "File $file exists."
        #     else
        #         echo "File $file does not exist."
        #         python -m cmass.survey.ngc_lightcone $postfix bias.hod.seed=$hod_seed survey.aug_seed=$aug_seed
        #     fi

        #     # diagnostics
        #     file=$outdir/$lhid/diag/lightcone/hod${hod_str}_aug${aug_str}.h5
        #     if [ -f $file ]; then
        #         echo "File $file exists."
        #     else
        #         echo "File $file does not exist."
        #         python -m cmass.diagnostics.summ_lc $postfix bias.hod.seed=$hod_seed survey.aug_seed=$aug_seed
        #     fi
        # done
    done
done