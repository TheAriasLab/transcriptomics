#!/bin/bash
#SBATCH -J load_sqlite_db_3-2
#SBATCH -o load_sqlite_db_3-2_out.log
#SBATCH -e load_sqlite_db_3-2_err.log
#SBATCH -n 32
#SBATCH --partition=normal

#cargar blast personalizado de cacao blastp

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_custom_blast --outfmt6 /BIOS-Share/home/areyesb/Anotacion/Blastp/blastp4.outfmt6 --prog blastp --dbtype blastp_cacao