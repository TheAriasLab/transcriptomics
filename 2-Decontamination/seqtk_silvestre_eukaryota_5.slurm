#!/bin/bash
#SBATCH -J seqtk_silvestre
#SBATCH -o seqtk_silvestre_out.log
#SBATCH -e seqtk_silvestre_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/seqtk/1.2

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217278_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217278_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217278_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217278_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217278_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217278_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217279_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217279_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217279_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217279_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217279_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217279_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217280_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217280_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217280_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217280_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217280_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217280_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217281_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217281_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217281_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217281_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217281_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217281_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217282_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217282_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217282_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217282_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217282_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217282_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217283_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217283_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217283_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217283_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217283_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217283_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217284_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217284_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217284_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217284_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217284_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217284_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217292_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217292_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217292_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217292_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217292_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217292_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217294_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217294_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217294_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217294_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217294_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217294_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217297_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217297_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217297_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217297_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217297_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre_filtrado_eukaryota/SRR3217297_2_Paired_filtrado.fastq