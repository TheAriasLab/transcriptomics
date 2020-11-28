#!/bin/bash
#SBATCH -J kaiju_cultivado
#SBATCH -o kaiju_cultivado.out
#SBATCH -e kaiju_cultivado.err
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate kaiju_env

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217276_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217276_2_Paired.fastq -o kaiju_cultivado1.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217277_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217277_2_Paired.fastq -o kaiju_cultivado2.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217298_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217298_2_Paired.fastq -o kaiju_cultivado3.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217299_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217299_2_Paired.fastq -o kaiju_cultivado4.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217301_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217301_2_Paired.fastq -o kaiju_cultivado5.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217304_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217304_2_Paired.fastq -o kaiju_cultivado6.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217315_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217315_2_Paired.fastq -o kaiju_cultivado7.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217317_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217317_2_Paired.fastq -o kaiju_cultivado8.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217318_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217318_2_Paired.fastq -o kaiju_cultivado9.out

kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217319_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217319_2_Paired.fastq -o kaiju_cultivado10.out