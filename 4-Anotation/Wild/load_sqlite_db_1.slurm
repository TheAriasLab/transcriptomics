#!/bin/bash
#SBATCH -J load_sqlite_db_1
#SBATCH -o load_sqlite_db_1_out.log
#SBATCH -e load_sqlite_db_1_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite init --gene_trans_map /BIOS-Share/home/areyesb/Anotacion/sqlite-db/Trinity.fasta.Silv.gene_trans_map --transcript_fasta /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta --transdecoder_pep /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder.pep