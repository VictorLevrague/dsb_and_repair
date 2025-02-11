#!/bin/sh

#SBATCH --job-name=DSB
#SBATCH --output="/sps/gdrmi2b/levrague/dsb_and_repair/log/%x_%j_%N.out"
#SBATCH --error="/sps/gdrmi2b/levrague/dsb_and_repair/log/%x.%j.err"
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8000
#SBATCH --time 7-00:00 # time (D-HH:MM)
#SBATCH --licenses=sps

export FILE=$1
export CopyNb=$2
export particle_energy=$3
export irradiation_id=$4
export dose=$5
export output_folder_name=output${particle_energy}_MeV/outputIrradiation${irradiation_id}/output_${dose}Gy

# Debugging statements
echo "Running test job with parameters:"
echo "  FILE: $FILE"
echo "  CopyNb: $CopyNb"
echo "  particle_energy: $particle_energy"
echo "  irradiation_id: $irradiation_id"
echo "  dose: $dose"
echo "  output_folder_name: $output_folder_name"

# Test script instead of dsbandrepair
output_file="/sps/gdrmi2b/levrague/dsb_and_repair/output/test_output_$SLURM_JOB_ID.txt"
echo "Test job started" > $output_file
echo "Processing file: $FILE" >> $output_file
sleep 5
echo "Test job completed successfully" >> $output_file

# Confirm completion
echo "Finished test job"
echo $FILE

