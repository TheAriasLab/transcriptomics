#!/bin/bash
#SBATCH -J report_trinotate
#SBATCH -o report_trinotate_out.log
#SBATCH -e report_trinotate_err.log
#SBATCH -n 32
#SBATCH --partition=normal

#generar reporte en xlsx

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite report [opts] > trinotate_annotation_report.xls