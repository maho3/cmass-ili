#!/bin/bash
#SBATCH --job-name=mtnglike_hm   # Job name
#SBATCH --time=02:00:00         # Time limit
#SBATCH --account=phy240043   # Account name
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=128            # Number of tasks
#SBATCH --partition=highmem      # Partition name

# fixing previous runs
# offset=0 401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,682
# offset=1000 430 457 481 548 605 662 682 733 734 735 736 737 738 739 740 741 742 743 744 745 746 747 748 749 750 751 752 753 754 755 756 757 759 760 942 969

# SLURM_ARRAY_TASK_ID=400
offset=0
task_ids=(401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421 422 423 424 425 426 427 428 429 430 682)
# task_ids=(430 457 481 548 605 662 682 733 734 735 736 737 738 739 740 741 742 743 744 745 746 747 748 749 750 751 752 753 754 755 756 757 759 760 942 969)


module restore cmass
conda activate cmassrun

# Command to run for each lhid
cd /home/x-mho1/git/ltu-cmass-run
outdir=/anvil/scratch/x-mho1/cmass-ili/mtnglike/fastpm/L3000-N384


for SLURM_ARRAY_TASK_ID in "${task_ids[@]}"; do
    echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID offset=$offset"

    lhid=$((SLURM_ARRAY_TASK_ID + offset))

    # check if nbody.h5 exists
    file=$outdir/$lhid/nbody.h5
    if [ -f $file ]; then
        echo "File $file exists."
    else
        echo "File $file does not exist."
        python -m cmass.nbody.fastpm nbody=mtnglike nbody.lhid=$lhid +nbody.postprocess=False

        sbatch --export=offset=$offset --array=$SLURM_ARRAY_TASK_ID jobs/slurm_mtng_postprocess.sh 
    fi
done