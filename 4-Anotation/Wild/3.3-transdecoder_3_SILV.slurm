#!/bin/bash
#SBATCH -J transdecoder_3_SILV
#SBATCH -o transdecoder_3_SILV_out.log
#SBATCH -e transdecoder_3_SILV_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/ncbi-blast/2.7.1

blastp -query /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep -db /BIOS-Share/home/areyesb/ProteCacao_Uniprot/uniprot-theobroma_cacao-filtered.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > blastpuniprotcacao.outfmt6