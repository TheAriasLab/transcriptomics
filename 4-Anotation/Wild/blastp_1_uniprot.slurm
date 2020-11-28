#!/bin/bash
#SBATCH -J blast_1_uniprot
#SBATCH -o blast_1_uniprot_out.log
#SBATCH -e blast_1_uniprot_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/ncbi-blast/2.7.1

blastx -query /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta -db /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx1.outfmt6