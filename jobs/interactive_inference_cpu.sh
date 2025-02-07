#!/bin/bash

module restore cmass_env
conda activate cmass

# exp_index=$SLURM_ARRAY_TASK_ID
exp_index=0

# Command to run for each lhid
cd /home/x-dbartlett/ltu-cmass

nbody=mtnglike
sim=fastpm
infer=bispectrum

halo=False
galaxy=False
ngc=False
sgc=False
mtng=True

#extras="nbody.zf=0.500015"
device=cpu

postfix="nbody=$nbody sim=$sim infer=$infer infer.exp_index=$exp_index"
postfix="$postfix infer.halo=$halo infer.galaxy=$galaxy"
postfix="$postfix infer.ngc_lightcone=$ngc infer.sgc_lightcone=$sgc infer.mtng_lightcone=$mtng"
postfix="$postfix infer.device=$device $extras"

echo "Running inference with $postfix"
python -m cmass.infer.preprocess $postfix
#python -m cmass.infer.train $postfix
