#!/bin/sh


### User parameters ###
export PATH=/sps/gdrmi2b/levrague/dsb_and_repair
export PARTICLE_ENERGY="H_150"
export IRRADIATION_ID=1 #1 or 2
export DOSE=1 #Gy
#######################

export split_folder_name=$PATH/split_phys/split${PARTICLE_ENERGY}_MeV/splitIrradiation${IRRADIATION_ID}/split_${DOSE}Gy

export job_prefix=$PARTICLE_ENERGY
export jobid="$job_prefix"
export job_path=$PATH/jobs
export jobname="$job_path"/"$job_prefix".sh
list_run=`ls $split_folder_name`

truncate -s 0 $jobname

for k in $list_run
do
	echo .$PATH/pbc $split_folder_name/$k 
	echo $split_folder_name/$k >> $jobname
	chmod u+x $jobname
done

export copy_nb=0

export log_dir=$PATH/log/log${PARTICLE_ENERGY}_MeV/splitIrradiation${IRRADIATION_ID}/split_${DOSE}Gy/split_phys

mkdir -p $log_dir
find $log_dir -type f -delete

for i in `cat $jobname`
    sbatch -J ${PARTICLE_ENERGY} --output=${log_dir}/${PARTICLE_ENERGY}_${copy_nb}.out $PATH/submit_slurm/1_phys/batch_parameters_phys_dsb_and_repair.sh $i $copy_nb $PARTICLE_ENERGY $IRRADIATION_ID $DOSE
	((copy_nb++))
done
