#!/bin/bash
#SBATCH -J blast_3_uniprot_refer
#SBATCH -o blast_3_uniprot_refer_out.log
#SBATCH -e blast_3_uniprot_refer_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/ncbi-blast/2.7.1

blastx -query /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta -db /BIOS-Share/home/areyesb/ProteCacao_Uniprot/uniprot-theobroma_cacao-filtered.fasta -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx3.outfmt6