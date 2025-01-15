#!/bin/sh
export path=/sps/gdrmi2b/levrague/dsb_and_repair

export particle_energy="C_290"
export irradiation_id=1 #1 or 2

export split_folder_name=split${particle_energy}_MeV/splitIrradiation${irradiation_id}

export job_prefix=$particle_energy
export jobid="$job_prefix"
export job_path=$path/jobs
export jobname="$job_path"/"$job_prefix".sh
list_run=`ls $path/split_phys/split${particle_energy}_MeV/splitIrradiation${irradiation_id}`

truncate -s 0 $jobname

for k in $list_run
do
	echo .$path/pbc $path/split_phys/split${particle_energy}_MeV/splitIrradiation${irradiation_id}/$k 
	echo $path/split_phys/split${particle_energy}_MeV/splitIrradiation${irradiation_id}/$k >> $jobname
	chmod u+x $jobname
done

export copy_nb=0

for i in `cat $jobname`
do
    mkdir -p $path/log/log${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_phys
    find $path/log/log${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_phys -type f -delete
	sbatch -J $path/log${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_phys/${job_prefix} batch_parameters_phys_dsb_and_repair.sh $i $copy_nb $particle_energy $irradiation_id
	((copy_nb++))
done
