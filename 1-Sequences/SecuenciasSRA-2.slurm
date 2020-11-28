#!/bin/bash
#SBATCH -J sec
#SBATCH -o sec.log
#SBATCH -e sec.log
#SBATCH -n 16
#SBATCH --partition=normal

module load software/bioinformatics/sratoolkit/2.8.1

#cambiar el tipo de archivo sra por uno que pueda leer trinity

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files *.sra -o /BIOS-Share/home/areyesb/Secuencias

 




