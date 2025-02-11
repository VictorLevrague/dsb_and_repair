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

export split_folder_name=$path/split_chem/split${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_${dose}Gy
export output_folder_name=$path/output/output${particle_energy}_MeV/outputIrradiation${irradiation_id}/output_${dose}Gy

mkdir -p $split_folder_name

for i in $(seq 0 $N)
do
  cp macros/run_chem.mac $split_folder_name/run_chem_$i.mac
done

for i in $(seq 0 $N)
do
  #L'emplacement de l'output est indiqué au moment de lancer le job : gros pb
  sed -i -e "s#/random/setSeeds 123456789 1#/random/setSeeds 123456789 $i#" $split_folder_name/run_chem_$i.mac
  mkdir -p $output_folder_name/Copy$i/chem_output
  echo $split_folder_name/run_chem_$i.mac
done
