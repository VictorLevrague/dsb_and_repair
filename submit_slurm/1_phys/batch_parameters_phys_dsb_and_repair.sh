#!/bin/sh

#SBATCH --job-name=DSB
#SBATCH --output="/sps/gdrmi2b/levrague/dsb_and_repair/log/%x-%j.out"
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 5000
#SBATCH --time 0-48:00 # time (D-HH:MM)
#SBATCH --licenses=sps

export FILE=$1
export CopyNb=$2

export particle_energy=$3
export irradiation_id=$4
export output_folder_name=output${particle_energy}_MeV/outputIrradiation${irradiation_id}

srun -n 1 --exclusive --cpus-per-task 1 /sps/gdrmi2b/levrague/dsb_and_repair/dsbandrepair $FILE phys $CopyNb /sps/gdrmi2b/levrague/dsb_and_repair/output/${output_folder_name} 1


echo $FILE
