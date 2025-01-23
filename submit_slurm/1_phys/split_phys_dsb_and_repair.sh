#!/bin/env bash


### User parameters ###
N=1000;
export PATH=/sps/gdrmi2b/levrague/dsb_and_repair
export PARTICLE_ENERGY="H_150"
export NB_PARTICLES=2976
export IRRADIATION_ID=1 #1 (HC nucleus) or 2 (nucleus with decompaction)
export DOSE=1 #Gy
#######################

N=`expr "$N" - 1`;

export split_folder_name=$PATH/split_phys/split${PARTICLE_ENERGY}_MeV/splitIrradiation${IRRADIATION_ID}/split_${DOSE}Gy
export output_folder_name=$PATH/output/output${PARTICLE_ENERGY}_MeV/outputIrradiation${IRRADIATION_ID}/output_${DOSE}Gy

mkdir -p $split_folder_name
mkdir -p $output_folder_name

for i in $(seq 0 $N)
do
  cp $PATH/macros/run_phys.mac $split_folder_name/run_phys_$i.mac
done

for i in $(seq 0 $N)
do
  #L'emplacement de l'output est indiqué au moment de lancer le job et pas dans la macro : gros pb
  echo $split_folder_name/run_phys_$i.mac
  sed -i -e "s#/random/setSeeds 123456789 1#/random/setSeeds 123456789 $i#" $split_folder_name/run_phys_$i.mac
  sed -i -e "s#/run/beamOn 1#/run/beamOn $NB_PARTICLES#" $split_folder_name/run_phys_$i.mac
  if [ $IRRADIATION_ID == 2 ]; then
    echo /sps/gdrmi2b/levrague/dsb_and_repair/geometry_creation/output/output${PARTICLE_ENERGY}_MeV/irradiated_geometry_file_${i}.fab2g4dna
    sed -i -e "s#/dsbandrepair/det/celldefinitionfile        /sps/gdrmi2b/levrague/dsb_and_repair/dnafabric_geometries/geometrie_fibroblast_hc_floriane.fab2g4dna#/dsbandrepair/det/celldefinitionfile        /sps/gdrmi2b/levrague/dsb_and_repair/geometry_creation/output/output${PARTICLE_ENERGY}_MeV/irradiated_geometry_file_${i}.fab2g4dna#" $split_folder_name/run_phys_$i.mac
  fi
  mkdir -p $output_folder_name/Copy$i
  mkdir -p $output_folder_name/Copy$i/phys_output
  mkdir -p $output_folder_name/Copy$i/chem_input
done
