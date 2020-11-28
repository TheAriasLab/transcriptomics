#!/bin/bash
#SBATCH -J blast_2_uniprot
#SBATCH -o blast_2_uniprot_out.log
#SBATCH -e blast_2_uniprot_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/ncbi-blast/2.7.1

blastp -query /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder.pep -db /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp2.outfmt6