#!/bin/bash
#SBATCH -J load_sqlite_6_signal
#SBATCH -o load_sqlite_6_signal_out.log
#SBATCH -e load_sqlite_6_signal_err.log
#SBATCH -n 32
#SBATCH --partition=normal

#cargar bases de peptidos señal por signalP

module load devtools/bioconda/bioconda3
source activate trinotate_env

Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_signalp /BIOS-Share/home/areyesb/Anotacion/signal/signalp.out