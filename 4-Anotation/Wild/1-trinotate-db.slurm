#!/bin/bash
#SBATCH -J trinotate_db
#SBATCH -o trinotate-db_out.log
#SBATCH -e trinotate-db_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate trinotate_env

Build_Trinotate_Boilerplate_SQLite_db.pl /BIOS-Share/home/areyesb/Anotacion/Trinotate-db

