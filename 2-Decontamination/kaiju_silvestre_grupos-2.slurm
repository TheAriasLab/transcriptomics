#!/bin/bash
#SBATCH -J kaiju_silvestre
#SBATCH -o kaiju_silvestre.out
#SBATCH -e kaiju_silvestre.err
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate kaiju_env

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217278_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217278_2_Paired.fastq -o kaiju_silvestre1.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217279_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217279_2_Paired.fastq -o kaiju_silvestre2.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217280_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217280_2_Paired.fastq -o kaiju_silvestre3.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217281_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217281_2_Paired.fastq -o kaiju_silvestre4.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217282_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217282_2_Paired.fastq -o kaiju_silvestre5.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217283_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217283_2_Paired.fastq -o kaiju_silvestre6.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217284_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217284_2_Paired.fastq -o kaiju_silvestre7.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217292_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217292_2_Paired.fastq -o kaiju_silvestre8.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217294_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217294_2_Paired.fastq -o kaiju_silvestre9.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217297_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217297_2_Paired.fastq -o kaiju_silvestre10.out