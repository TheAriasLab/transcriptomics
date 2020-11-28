#!/bin/bash
#SBATCH -J seqtk_cultivado_eu
#SBATCH -o seqtk_cultivado_eu_out.log
#SBATCH -e seqtk_cultivado_eu_err.log
#SBATCH -n 32
#SBATCH --partition=normal

module load software/bioinformatics/seqtk/1.2

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217276_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217276_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217276_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217276_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217276_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217276_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217277_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217277_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217277_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217277_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217277_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217277_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217298_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217298_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217298_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217298_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217298_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217298_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217299_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217299_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217299_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217299_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217299_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217299_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217301_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217301_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217301_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217301_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217301_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217301_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217304_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217304_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217304_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217304_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217304_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217304_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217315_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217315_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217315_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217315_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217315_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217315_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217317_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217317_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217317_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217317_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217317_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217317_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217318_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217318_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217318_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217318_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217318_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217318_2_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217319_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217319_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217319_1_Paired_filtrado.fastq

seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217319_2_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_cult-Eukaryota_SRR3217319_2.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado_eukaryota/SRR3217319_2_Paired_filtrado.fastq