#!/bin/bash
#SBATCH -J transdecoder_2_SILV
#SBATCH -o transdecoder_2_SILV_out.log
#SBATCH -e transdecoder_2_SILV_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate transdecoder_env

TransDecoder.Predict -t /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta