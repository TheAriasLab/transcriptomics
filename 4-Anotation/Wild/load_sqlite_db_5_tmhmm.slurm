#!/bin/bash
#SBATCH -J load_sqlite_5_tmhmm
#SBATCH -o load_sqlite_5_tmhmm_out.log
#SBATCH -e load_sqlite_5_tmhmm_err.log
#SBATCH -n 32
#SBATCH --partition=normal

#cargar bases de dominio transmembranas generados por tmhmm

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_tmhmm /BIOS-Share/home/areyesb/Anotacion/tmhmm/tmhmm.out