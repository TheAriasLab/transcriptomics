#!/bin/bash
#SBATCH -J load_sqlite_db_2-2
#SBATCH -o load_sqlite_db_2-2_out.log
#SBATCH -e load_sqlite_db_2-2_err.log
#SBATCH -n 32
#SBATCH --partition=normal

#cargar blast no personalizados blastx

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_swissprot_blastx /BIOS-Share/home/areyesb/Anotacion/Blastp/blastx1.outfmt6