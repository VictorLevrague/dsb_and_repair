#!/bin/sh

export PATH=/sps/gdrmi2b/levrague/dsb_and_repair
export JOB_PREFIX="DSB"
export PARTICLE_ENERGY="O_350"
export IRRADIATION_ID=2 #1 or 2
export DOSE=1 #Gy
N=10;

for i in $(seq 0 $N)
do
    sbatch -J log${PARTICLE_ENERGY}_MeV/splitIrradiation${IRRADIATION_ID}/split_${DOSE}Gy/split_ana/${JOB_PREFIX} $PATH/submit_slurm/3_analysis/batch_parameters_analysis_dsb_and_repair.sh $PATH/macros/analysis.in $i $PARTICLE_ENERGY $IRRADIATION_ID $DOSE
done
