#!/bin/bash
#SBATCH -J transdecoder_3.1_SILV
#SBATCH -o transdecoder_3.1_SILV_out.log
#SBATCH -e transdecoder_3.1_SILV_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate hmmer_env

hmmpress /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Pfam-A.hmm