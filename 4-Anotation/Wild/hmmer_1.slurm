#!/bin/bash
#SBATCH -J hmmer_1
#SBATCH -o hmmer_1_out.log
#SBATCH -e hmmer_1_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate hmmer_env

hmmscan --cpu 32 --domtblout /BIOS-Share/home/areyesb/Anotacion/Transdecoder/pfam.domtblout /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Pfam-A.hmm /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder.pep > pfam.log 	