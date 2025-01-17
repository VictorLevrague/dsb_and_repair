#!/bin/env bash

##############################################################
# This script split "run.mac" (example for 10 jobs). Execute with:
# ./split_mac_file.sh

export path=/sps/gdrmi2b/levrague/dsb_and_repair

N=1000;
export particle_energy="H_150"
export nb_particles=2976
export irradiation_id=1 #1 (HC nucleus) or 2 (nucleus with decompaction)
export dose=1 #Gy

N=`expr "$N" - 1`;

export split_folder_name=$path/split_phys/split${particle_energy}_MeV/splitIrradiation${irradiation_id}/split_${dose}Gy
export output_folder_name=$path/output/output${particle_energy}_MeV/outputIrradiation${irradiation_id}/output_${dose}Gy

mkdir -p $split_folder_name
mkdir -p $output_folder_name

for i in $(seq 0 $N)
do
  cp $path/macros/run_phys.mac $split_folder_name/run_phys_$i.mac
done

for i in $(seq 0 $N)
do
  #L'emplacement de l'output est indiqué au moment de lancer le job et pas dans la macro : gros pb
  echo $split_folder_name/run_phys_$i.mac
  sed -i -e "s#/random/setSeeds 123456789 1#/random/setSeeds 123456789 $i#" $split_folder_name/run_phys_$i.mac
  sed -i -e "s#/run/beamOn 1#/run/beamOn $nb_particles#" $split_folder_name/run_phys_$i.mac
  if [ $irradiation_id == 2 ]; then
    echo /sps/gdrmi2b/levrague/dsb_and_repair/geometry_creation/output/output${particle_energy}_MeV/irradiated_geometry_file_${i}.fab2g4dna
    sed -i -e "s#/dsbandrepair/det/celldefinitionfile        /sps/gdrmi2b/levrague/dsb_and_repair/dnafabric_geometries/geometrie_fibroblast_hc_floriane.fab2g4dna#/dsbandrepair/det/celldefinitionfile        /sps/gdrmi2b/levrague/dsb_and_repair/geometry_creation/output/output${particle_energy}_MeV/irradiated_geometry_file_${i}.fab2g4dna#" $split_folder_name/run_phys_$i.mac
  fi
  mkdir -p $output_folder_name/Copy$i
  mkdir -p $output_folder_name/Copy$i/phys_output
  mkdir -p $output_folder_name/Copy$i/chem_input
done
