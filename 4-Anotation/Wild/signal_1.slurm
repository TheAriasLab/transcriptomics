#!/bin/bash
#SBATCH -J signal
#SBATCH -o signal_1_out.log
#SBATCH -e signal_1_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/signalp/4.1

signalp -f short -n signalp.out /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Trinity.fasta.transdecoder.pep