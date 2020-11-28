#!/bin/bash
#SBATCH -J trim
#SBATCH -o trim_out.log
#SBATCH -e trim_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/trimmomatic/0.36

#cultivado
#SRR3217276 SRR3217277 SRR3217298 SRR3217299 SRR3217301 SRR3217304 SRR3217315 SRR3217317 SRR3217318 SRR3217319


TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217276_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217276_2.fastq SRR3217276_1_Paired.fastq SRR3217276_1_UNPaired.fastq SRR3217276_2_Paired.fastq SRR3217276_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217277_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217277_2.fastq SRR3217277_1_Paired.fastq SRR3217277_1_UNPaired.fastq SRR3217277_2_Paired.fastq SRR3217277_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217298_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217298_2.fastq SRR3217298_1_Paired.fastq SRR3217298_1_UNPaired.fastq SRR3217298_2_Paired.fastq SRR3217298_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217299_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217299_2.fastq SRR3217299_1_Paired.fastq SRR3217299_1_UNPaired.fastq SRR3217299_2_Paired.fastq SRR3217299_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217301_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217301_2.fastq SRR3217301_1_Paired.fastq SRR3217301_1_UNPaired.fastq SRR3217301_2_Paired.fastq SRR3217301_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217304_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217304_2.fastq SRR3217304_1_Paired.fastq SRR3217304_1_UNPaired.fastq SRR3217304_2_Paired.fastq SRR3217304_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217315_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217315_2.fastq SRR3217315_1_Paired.fastq SRR3217315_1_UNPaired.fastq SRR3217315_2_Paired.fastq SRR3217315_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217317_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217317_2.fastq SRR3217317_1_Paired.fastq SRR3217317_1_UNPaired.fastq SRR3217317_2_Paired.fastq SRR3217317_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217318_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217318_2.fastq SRR3217318_1_Paired.fastq SRR3217318_1_UNPaired.fastq SRR3217318_2_Paired.fastq SRR3217318_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217319_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217319_2.fastq SRR3217319_1_Paired.fastq SRR3217319_1_UNPaired.fastq SRR3217319_2_Paired.fastq SRR3217319_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3


