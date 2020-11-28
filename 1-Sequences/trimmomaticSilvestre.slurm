#!/bin/bash
#SBATCH -J trim
#SBATCH -o trim_out.log
#SBATCH -e trim_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/trimmomatic/0.36

#silvestre
#SRR3217278 SRR3217279 SRR3217280 SRR3217281 SRR3217282 SRR3217283 SRR3217284 SRR3217292 SRR3217294 SRR3217297

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217278_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217278_2.fastq SRR3217278_1_Paired.fastq SRR3217278_1_UNPaired.fastq SRR3217278_2_Paired.fastq SRR3217278_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217279_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217279_2.fastq SRR3217279_1_Paired.fastq SRR3217279_1_UNPaired.fastq SRR3217279_2_Paired.fastq SRR3217279_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217280_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217280_2.fastq SRR3217280_1_Paired.fastq SRR3217280_1_UNPaired.fastq SRR3217280_2_Paired.fastq SRR3217280_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217281_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217281_2.fastq SRR3217281_1_Paired.fastq SRR3217281_1_UNPaired.fastq SRR3217281_2_Paired.fastq SRR3217281_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217282_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217282_2.fastq SRR3217282_1_Paired.fastq SRR3217282_1_UNPaired.fastq SRR3217282_2_Paired.fastq SRR3217282_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217283_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217283_2.fastq SRR3217283_1_Paired.fastq SRR3217283_1_UNPaired.fastq SRR3217283_2_Paired.fastq SRR3217283_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217284_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217284_2.fastq SRR3217284_1_Paired.fastq SRR3217284_1_UNPaired.fastq SRR3217284_2_Paired.fastq SRR3217284_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217292_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217292_2.fastq SRR3217292_1_Paired.fastq SRR3217292_1_UNPaired.fastq SRR3217292_2_Paired.fastq SRR3217292_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217294_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217294_2.fastq SRR3217294_1_Paired.fastq SRR3217294_1_UNPaired.fastq SRR3217294_2_Paired.fastq SRR3217294_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3

TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217297_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217297_2.fastq SRR3217297_1_Paired.fastq SRR3217297_1_UNPaired.fastq SRR3217297_2_Paired.fastq SRR3217297_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3


