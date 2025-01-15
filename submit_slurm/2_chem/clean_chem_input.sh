#!/bin/env bash

##############################################################
# This script split "run.mac" (example for 10 jobs). Execute with:
# ./split_mac_file.sh

export path=/sps/gdrmi2b/levrague/dsb_and_repair

N=1000;
export particle_energy="O_350"

N=`expr "$N" - 1`;

for i in $(seq 0 $N)
do
  export output_folder_name=$path/output/output${particle_energy}_MeV/Copy$i/chem_input
  echo $output_folder_name
  list_run=`ls $output_folder_name`
  rm -rf $output_folder_name
done
