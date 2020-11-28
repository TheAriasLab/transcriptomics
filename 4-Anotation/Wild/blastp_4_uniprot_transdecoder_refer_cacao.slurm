#!/bin/bash
#SBATCH -J blast_4_uniprot_trans_refer
#SBATCH -o blast_4_uniprot_trans_refer_out.log
#SBATCH -e blast_4_uniprot_trans_refer_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/ncbi-blast/2.7.1

blastp -query /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder.pep -db /BIOS-Share/home/areyesb/ProteCacao_Uniprot/uniprot-theobroma_cacao-filtered.fasta -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp4.outfmt6