#!/bin/env bash

##############################################################
# This script split "run.mac" (example for 10 jobs). Execute with:
# ./split_mac_file.sh

N=1000;
export particle_energy="C_290"
export nb_particles=125
export irradiation_id=1 #1 (HC) or 2 (decompaction)

N=`expr "$N" - 1`;

export split_folder_name=split${particle_energy}_MeV/splitIrradiation${irradiation_id}
export output_folder_name=output${particle_energy}_MeV/outputIrradiation${irradiation_id}

mkdir -p split_phys/$split_folder_name
mkdir -p output/$output_folder_name

for i in $(seq 0 $N)
do
  cp macros/run_phys.mac split_phys/$split_folder_name/run_phys_$i.mac
done

for i in $(seq 0 $N)
do
  #L'emplacement de l'output est indiqué au moment de lancer le job : gros pb
  echo split_phys/$split_folder_name/run_phys_$i.mac
  sed -i -e "s#/random/setSeeds 123456789 1#/random/setSeeds 123456789 $i#" split_phys/$split_folder_name/run_phys_$i.mac
  sed -i -e "s#/run/beamOn 1#/run/beamOn $nb_particles#" split_phys/$split_folder_name/run_phys_$i.mac
  if [ $irradiation_id == 2 ]; then
    echo /sps/gdrmi2b/levrague/dsb_and_repair/geometry_creation/output/output${particle_energy}_MeV/irradiated_geometry_file_${i}.fab2g4dna
    sed -i -e "s#/dsbandrepair/det/celldefinitionfile        /sps/gdrmi2b/levrague/dsb_and_repair/dnafabric_geometries/geometrie_fibroblast_hc_floriane.fab2g4dna#/dsbandrepair/det/celldefinitionfile        /sps/gdrmi2b/levrague/dsb_and_repair/geometry_creation/output/output${particle_energy}_MeV/irradiated_geometry_file_${i}.fab2g4dna#" split_phys/$split_folder_name/run_phys_$i.mac
  fi
  mkdir -p output/$output_folder_name/Copy$i
  mkdir -p output/$output_folder_name/Copy$i/phys_output
  mkdir -p output/$output_folder_name/Copy$i/chem_input
done
