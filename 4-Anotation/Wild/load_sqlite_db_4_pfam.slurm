#!/bin/bash
#SBATCH -J load_sqlite_4_pfam
#SBATCH -o load_sqlite_4_pfam_out.log
#SBATCH -e load_sqlite_4_pfam_err.log
#SBATCH -n 32
#SBATCH --partition=normal

#cargar bdominio generados por pfam

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_pfam /BIOS-Share/home/areyesb/Anotacion/Transdecoder/pfam.domtblout