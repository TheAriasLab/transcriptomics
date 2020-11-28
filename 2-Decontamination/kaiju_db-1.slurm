#!/bin/bash

#SBATCH -J kaiju_db
#SBATCH -o kaiju_db.out
#SBATCH -e kaiju_db.err
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate kaiju_env

kaiju-makedb -s nr_euk