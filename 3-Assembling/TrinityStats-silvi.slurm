#!/bin/bash
#SBATCH -J TrinityStats-silvi
#SBATCH -o TrinityStats-silvi.log
#SBATCH -e TrinityStats-silvi.err
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/trinityrnaseq/2.8.4

TrinityStats.pl /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta