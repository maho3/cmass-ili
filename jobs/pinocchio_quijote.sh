#!/bin/sh -l
# FILENAME:  pinocchio_quijote_0

#SBATCH -A phy240043
#SBATCH -p shared # the default queue is "shared" queue
#SBATCH --nodes=1
#SBATCH --ntasks=32 
#SBATCH --time=1:00:00
#SBATCH --job-name pinocchio
#SBATCH --output=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out  # Output file for each array task
#SBATCH --error=/anvil/scratch/x-mho1/jobout/%x_%A_%a.out   # Error file for each array task

module purge 
module restore pinocchio
module list 
conda activate cmassrun

# Print the hostname of the compute node on which this job is running.
hostname

# Change to correct directory
cd /home/x-mho1/git/ltu-cmass-run
pwd

# Define the range for lhid
start_lhid=0
end_lhid=1

# Loop over lhid from start_lhid to end_lhid
for lhid in $(seq $start_lhid $end_lhid)
do
  python -m cmass.nbody.pinocchio nbody=pinocchio nbody.lhid=$lhid
done

exit 0
