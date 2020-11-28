#!/bin/bash
#SBATCH -J transdecoder_4_SILV
#SBATCH -o transdecoder_4_SILV_out.log
#SBATCH -e transdecoder_4_SILV_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate transdecoder_env

TransDecoder.Predict -t /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta --retain_pfam_hits /BIOS-Share/home/areyesb/Anotacion/Transdecoder/pfam.domtblout --retain_blastp_hits /BIOS-Share/home/areyesb/Anotacion/Transdecoder/blastpuniprotcacao.outfmt6