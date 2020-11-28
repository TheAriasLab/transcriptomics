#!/bin/bash
#SBATCH -J trinotate_db2
#SBATCH -o trinotate-db2_out.log
#SBATCH -e trinotate-db2_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate trinotate_env

makeblastdb -in uniprot_sprot.pep -dbtype prot