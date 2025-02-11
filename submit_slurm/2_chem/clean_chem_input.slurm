#!/bin/bash
#SBATCH -o logClean.txt
#SBATCH -e logCleanError.txt
#SBATCH -t 1-00:00
#SBATCH -n 1
#SBATCH --mem=8G

submit_slurm/2_chem/clean_chem_input.sh
