#!/bin/bash
#SBATCH -J tmhmm_1
#SBATCH -o tmhmm_1_out.log
#SBATCH -e tmhmm_1_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/tmhmm/2.0c

tmhmm --short < /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder.pep > tmhmm.out

	