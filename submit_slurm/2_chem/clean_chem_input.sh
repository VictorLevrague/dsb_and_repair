#!/bin/env bash

##############################################################
# This script split "run.mac" (example for 10 jobs). Execute with:
# ./split_mac_file.sh

export path=/sps/gdrmi2b/levrague/dsb_and_repair

N=1000;
export particle_energy="C_290"
export irradiation_id=1 #1 or 2
export dose=1 #Gy

N=`expr "$N" - 1`;

for i in $(seq 0 $N)
do
  export output_folder_name=$path/output/output${particle_energy}_MeV/outputIrradiation${irradiation_id}/output_${dose}Gy/Copy$i/chem_input
  echo $output_folder_name
  ls "$output_folder_name" 2>/dev/null
  rm -rf "$output_folder_name"
done
