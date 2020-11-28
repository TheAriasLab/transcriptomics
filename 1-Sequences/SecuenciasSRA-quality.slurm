#!/bin/bash
#SBATCH -J sec-quality1
#SBATCH -o sec-quality1.log
#SBATCH -e sec-quality1.log
#SBATCH -n 16
#SBATCH --partition=normal

module load software/bioinformatics/fastqc/0.11.4

#Correr para todas las secuencias
fastqc *.fastq






