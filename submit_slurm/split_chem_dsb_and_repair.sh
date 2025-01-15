#!/bin/env bash

##############################################################
# This script split "run.mac" (example for 10 jobs). Execute with:
# ./split_mac_file.sh

N=1;
export particle_energy="O_350"
export irradiation_id=2 #1 or 2

N=`expr "$N" - 1`;

export split_folder_name=split${particle_energy}_MeV/splitIrradiation${irradiation_id}
export output_folder_name=output${particle_energy}_MeV/outputIrradiation${irradiation_id}

mkdir -p split_chem/$split_folder_name

for i in $(seq 0 $N)
do
  cp macros/run_chem.mac split_chem/$split_folder_name/run_chem_$i.mac
done

for i in $(seq 0 $N)
do
  #L'emplacement de l'output est indiqué au moment de lancer le job : gros pb
  sed -i -e "s#/random/setSeeds 123456789 1#/random/setSeeds 123456789 $i#" split_chem/$split_folder_name/run_chem_$i.mac
  mkdir -p output/$output_folder_name/Copy$i/chem_output
done
