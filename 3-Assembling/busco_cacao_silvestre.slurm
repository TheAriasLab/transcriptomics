#!/bin/bash

#SBATCH -J busco_silv_euk
#SBATCH -o busco_silv_euk.out
#SBATCH -e busco_silv_euk.err
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate busco_env

busco -i /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta -o busco.out -l embryophyta_odb10 -m tran -c 32 -f