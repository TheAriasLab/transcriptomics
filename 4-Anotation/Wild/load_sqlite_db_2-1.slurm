#!/bin/bash
#SBATCH -J load_sqlite_db_2-1
#SBATCH -o load_sqlite_db_2-1_out.log
#SBATCH -e load_sqlite_db_2-1_err.log
#SBATCH -n 32
#SBATCH --partition=normal

#cargar blast no personalizados blastp

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_swissprot_blastp /BIOS-Share/home/areyesb/Anotacion/Blastp/blastp2.outfmt6