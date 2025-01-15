#!/bin/sh
export path=/sps/gdrmi2b/levrague/dsb_and_repair
export job_prefix="DSB"

export particle_energy="O_350"
export irradiation_id=2 #1 or 2
N=10;

for i in $(seq 0 $N)
do
    sbatch -J log${particle_energy}_MeV/splitIrradiation${irradiation_id}/${job_prefix} batch_parameters_analysis_dsb_and_repair.sh macros/analysis.in $i $particle_energy $irradiation_id
done
