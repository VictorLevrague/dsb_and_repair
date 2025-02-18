#!/bin/sh
export path=/sps/gdrmi2b/levrague/dsb_and_repair

export particle_energy="C_290"
export irradiation_id=1 #1 or 2
export dose=1 #Gy

export split_folder_name=$path/split_chem/split${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_${dose}Gy

export job_prefix=$particle_energy
export jobid="$job_prefix"
export job_path=$path/jobs
export jobname="$job_path"/"$job_prefix".sh
list_run=`ls $split_folder_name`

truncate -s 0 $jobname

for k in $list_run
do
	echo .$path/pbc $split_folder_name/$k 
	echo $split_folder_name/$k >> $jobname
	chmod u+x $jobname
done

export copy_nb=0

export log_dir=$path/log/log${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_${dose}Gy/split_chem

mkdir -p $log_dir
find $log_dir -type f -delete

for i in `cat $jobname`
do
	sbatch -J ${particle_energy} --output=${log_dir}/${particle_energy}_${copy_nb}.out $path/submit_slurm/2_chem/batch_parameters_chem_dsb_and_repair.sh $i $copy_nb $particle_energy $irradiation_id $dose
	((copy_nb++))
done
