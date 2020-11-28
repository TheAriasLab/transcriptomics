#!/bin/bash
#SBATCH -J transdecoder_3.2_SILV
#SBATCH -o transdecoder_3.2_SILV_out.log
#SBATCH -e transdecoder_3.2_SILV_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate hmmer_env

hmmscan --cpu 32 --domtblout pfam.domtblout /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Pfam-A.hmm /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep