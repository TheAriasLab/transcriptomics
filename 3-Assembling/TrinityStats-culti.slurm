#!/bin/bash
#SBATCH -J TrinityStats-culti
#SBATCH -o TrinityStats-culti.log
#SBATCH -e TrinityStats-culti.err
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/trinityrnaseq/2.8.4

TrinityStats.pl /BIOS-Share/home/areyesb/Trinity/Cultivado-Ensam-eukaryota/trinity_out_dir/Trinity.fasta
