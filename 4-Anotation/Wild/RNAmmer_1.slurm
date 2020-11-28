#!/bin/bash
#SBATCH -J RNAmmer
#SBATCH -o RNAmmer_1_out.log
#SBATCH -e RNAmmer_1_err.log
#SBATCH -n 32
#SBATCH --partition=normal



module load devtools/bioconda/bioconda3
source activate trinotate_env

/BIOS-Share/Software/devtools/bioconda/miniconda3/python3.7/envs/trinotate_env/bin/RnammerTranscriptome.pl --transcriptome /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam-eukaryota/trinity_out_dir/Trinity.fasta --path_to_rnammer /BIOS-Share/home/areyesb/Anotacion/RNAmmer/rnammer/rnammer