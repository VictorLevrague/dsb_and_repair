#!/bin/sh
export path=/sps/gdrmi2b/levrague/dsb_and_repair
export job_prefix="DSB"

export particle_energy="O_350"
export irradiation_id=2 #1 or 2

export jobid="$job_prefix"
export job_path=$path/jobs
export jobname="$job_path"/"$job_prefix".sh
list_run=`ls $path/split_chem/split${particle_energy}_MeV/splitIrradiation${irradiation_id}`

truncate -s 0 $jobname

for k in $list_run
do
	echo .$path/pbc $path/split_chem/split${particle_energy}_MeV/splitIrradiation${irradiation_id}/$k 
	echo $path/split_chem/split${particle_energy}_MeV/splitIrradiation${irradiation_id}/$k >> $jobname
	chmod u+x $jobname
done

export copy_nb=0

for i in `cat $jobname`
do
    mkdir -p $path/log/log${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_chem
    find $path/log/log${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_chem -type f -delete
	sbatch -J $path/log${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_chem/${job_prefix} batch_parameters_chem_dsb_and_repair.sh $i $copy_nb $particle_energy $irradiation_id
	((copy_nb++))
done
